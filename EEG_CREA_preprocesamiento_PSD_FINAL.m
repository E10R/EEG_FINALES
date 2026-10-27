%% EEG_CREA - PREPROCESAMIENTO v5 + ECG-QRS QC | 26-oct-2026
% Pipeline: EDF -> excluir ECG -> 0.5-40 Hz -> detectar canales malos
% SIN ASR -> quitar canales malos -> UN ASR con rechazo -> ICA+ICLabel
% -> interpolacion -> guardar .set/.mat y trazabilidad.
%
% IMPORTANTE: v3 candidato a validar en pilotos ANTES de analizar toda la muestra.
% No calcula PSD / FOOOF / PLV ni decide inclusion de ventanas.
% No modifica EDF originales. Requiere EEGLAB, clean_rawdata, firfilt,
% ICLabel y DIPFIT (para plantilla estandar 10-20).
% La tasa de muestreo EDF se mantiene (esperada: 512 Hz).
%
% Dos fases en 1 y 2 se descubren recorriendo subcarpetas y *.edf;
% no se presuponen nombres de fases ni numero exacto de EDF.
clearvars; clc;

%% CONFIGURACION
EEGLAB_DIR = 'C:\Users\evadi\OneDrive\Escritorio\eeglab_current\eeglab2026.0.0';
ROOT = 'C:\Users\evadi\OneDrive\Escritorio\FENNSI\CREATIVIDAD_DEFINITIVOS_EEG\Archivos_EEG';
DEST = fullfile(ROOT,'EEG_preprocesados'); % nombre pedido, conservando grafia literal
GROUPS = {'1','2'};
% {} = todos los sujetos; para prueba, por ejemplo: {'S07'}
SUBJECT_FILTER = {'S17'}; % {} = TODOS los sujetos
PHASE_FILTER = {'tSMS'}; % {} = TODAS las fases; {'tSMS'} = solo tSMS
% false = NO sobrescribir archivos existentes. Recomendado.
OVERWRITE = false;
SAVE_FULL_MAT = true;    % ademas del .set, guarda EEG_preprocess_*.mat
SAVE_FILTERED_QC = true; % true: .set adicional filtrado antes de ASR (pesado)
SAVE_PRE_ICA_QC = true;  % true: .set adicional despues ASR y antes de ICA
SAVE_ECG_REFERENCE = true;
ECG_QRS_QC = true; % deteccion de R e ICs sincronizadas, SOLO descriptiva, sin rechazo automatico
SAVE_PSD_QC = true; % Figura ANTES vs DESPUES por sujeto/fase
VERIFY_SAVED_FILES = true; % REABRIR .set, .mat y mapa para comparar exactamente
REQUIRE_NO_F3 = true; % EEG_CREA: 19 electrodos EEG, F3 ausente
PSD_QC_WINDOW_SEC = 2; % Welch SOLO en fragmentos continuos originales
PSD_QC_YLIM_COMMON = true; % Eje Y compartido en los dos paneles


% PARAMETROS ORIGINALES
FS_EXPECTED = 512;
FILT_LOW = 0.5;
FILT_HIGH = 40;
FLATLINE = 5;
CHAN_CORR = 0.8;
LINE_NOISE = 4;
BURST_CRITERION = 20;
WINDOW_CRITERION = 0.25;
WINDOW_TOLERANCES = [-Inf 7];
ASR_DISTANCE = 'Euclidian'; % grafia tal como en codigo original
ICLABEL_THRESHOLD = 0.70; % muscle eye heart line channel; NO Other
ASR_MAXMEM = 64; % igual que en el piloto; documentado para repetibilidad
QC_EPOCH_SECONDS = 120; % informacion descriptiva; NO segmenta el EEG
MIN_WELCH_WINDOWS_2MIN = 15; % provisional: 15 ventanas completas de 2 s (30 s efectivos)
MIN_GOOD_CHANNELS = 8;

%% INICIALIZACION Y VERIFICACIONES
assert(isfile(fullfile(EEGLAB_DIR,'eeglab.m')),'No se encontro eeglab.m: %s',EEGLAB_DIR);
assert(isfolder(ROOT),'Carpeta raiz inexistente: %s',ROOT);
addpath(EEGLAB_DIR);
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab('nogui'); %#ok<ASGLU>
required = {'pop_biosig','pop_chanedit','pop_eegfiltnew','pop_clean_rawdata', ...
    'pop_runica','pop_iclabel','pop_subcomp','pop_interp','pop_saveset'};
for i=1:numel(required)
    assert(~isempty(which(required{i})),'No esta disponible: %s',required{i});
end
LOC = fullfile(fileparts(which('eeglab')),'plugins','dipfit', ...
    'standard_BEM','elec','standard_1020.elc');
assert(isfile(LOC),'No encuentro standard_1020.elc en %s',LOC);
if ~isfolder(DEST), mkdir(DEST); end
runstamp = datestr(now,'yyyymmdd_HHMMSS');
QCDIR = fullfile(DEST,['AUDITORIA_' runstamp]);
mkdir(QCDIR);

% Log de texto (diary puede cerrarse al acabar)
diary(fullfile(QCDIR,'CONSOLa_EJECUCION.txt'));
cleanupDiary = onCleanup(@()diary('off')); %#ok<NASGU>
fprintf('Inicio: %s\nOrigen: %s\nDestino: %s\nEEGLAB: %s\n',datestr(now),ROOT,DEST,which('eeglab'));
fprintf('ASR: BC=%g, WC=%g, rechazo=on | ICLabel %.2f\n', ...
    BURST_CRITERION,WINDOW_CRITERION,ICLABEL_THRESHOLD);
fprintf('WELCH_2MIN: ventanas de %.0f s, minimo provisional %d ventanas; sin exclusion global.\n', ...
    PSD_QC_WINDOW_SEC,MIN_WELCH_WINDOWS_2MIN);

SUMMARY = table(); CHANNELS = table(); COMPONENTS = table(); QRS_SUMMARY=table();
WINDOWS = table(); FRAGMENTS = table(); INCIDENTS = table();

%% RECORRER GRUPOS / SUJETOS / EDF
for gi=1:numel(GROUPS)
    grp = GROUPS{gi};
    idir = fullfile(ROOT,grp);
    if ~isfolder(idir)
        INCIDENTS = [INCIDENTS; table(string(grp),"","", "Grupo inexistente", ...
            'VariableNames',{'Grupo','Sujeto','Fase','Incidencia'})]; %#ok<AGROW>
        continue;
    end
    dsub = dir(idir);
    dsub = dsub([dsub.isdir] & ~ismember({dsub.name},{'.','..'}));
    for si=1:numel(dsub)
        subject = dsub(si).name;
        if ~isempty(SUBJECT_FILTER) && ~any(strcmpi(subject,SUBJECT_FILTER)), continue; end
        srcSub = fullfile(idir,subject);
        EDFs = dir(fullfile(srcSub,'*.edf'));
        % Soporta EDF con extension mayuscula (evita repeticion en Windows)
        if isempty(EDFs), EDFs = dir(fullfile(srcSub,'*.EDF')); end
        for fi=1:numel(EDFs)
            [~,stem] = fileparts(EDFs(fi).name);
            if ~isempty(PHASE_FILTER) && ~any(contains(stem,PHASE_FILTER,'IgnoreCase',true))
                continue;
            end
            srcEDF = fullfile(srcSub,EDFs(fi).name);
            outSub = fullfile(DEST,grp,subject);
            if ~isfolder(outSub), mkdir(outSub); end
            outfile = fullfile(outSub,[stem '_preprocesado.set']);
            fprintf('\n===== %s | %s | %s =====\n',grp,subject,stem);
            try
                completeMarker = fullfile(outSub,[stem '_QC_COMPLETO.ok']);
                if isfile(outfile) && ~OVERWRITE
                    if ~isfile(completeMarker)
                        error(['Hay un .set existente sin marcador de QC completo: %s. ' ...
                            'No se sobrescribe automaticamente; revisar o mover el dataset.'],outfile);
                    end
                    fprintf('OMITIDO (QC completo): existe %s\n',outfile);
                    INCIDENTS = [INCIDENTS; table(string(grp),string(subject),string(stem), ...
                        "Archivo existente; no sobrescrito",'VariableNames', ...
                        {'Grupo','Sujeto','Fase','Incidencia'})]; %#ok<AGROW>
                    continue;
                end

                %% 1) IMPORTAR EDF Y SEPARAR ECG
                EEG = pop_biosig(srcEDF);
                EEG = eeg_checkset(EEG);
                assert(EEG.trials==1,'Se esperaba EEG continuo; trials=%d',EEG.trials);
                assert(abs(EEG.srate-FS_EXPECTED)<1e-6, ...
                    'Fs inesperada: %.3f Hz',EEG.srate);
                fs=EEG.srate;
                rawlabels=string({EEG.chanlocs.labels});
                ecgIdx=find(strcmpi(rawlabels,'ECG1+') | startsWith(upper(rawlabels),'ECG'));
                auxIdx=find(contains(upper(rawlabels),'MKR'));
                assert(numel(ecgIdx)<=1,'Se encontraron varios canales ECG; revisar: %s',stem);
                % ECG no entra en EEG/ASR/ICA. Se conserva separado para QC de ICs.
                ecg_orig = [];
                if numel(ecgIdx)==1, ecg_orig=double(EEG.data(ecgIdx,:)); end
                auxiliaryIdx=unique([ecgIdx(:);auxIdx(:)]);
                if ~isempty(auxiliaryIdx)
                    EEG=pop_select(EEG,'nochannel',cellstr(rawlabels(auxiliaryIdx)));
                end
                assert(EEG.nbchan==19,'Esperaba 19 canales EEG tras quitar auxiliares, hay %d',EEG.nbchan);
                if REQUIRE_NO_F3
                    assert(~any(strcmpi(string({EEG.chanlocs.labels}),'F3')), ...
                        'Se encontro F3: el montaje EEG_CREA debe tener 19 EEG sin F3.');
                end
                EEG=pop_chanedit(EEG,'lookup',LOC);
                EEG=eeg_checkset(EEG);
                referenceChanlocs=EEG.chanlocs; % LOCS antes de eliminar canales
                originalLabels=string({EEG.chanlocs.labels});
                N=EEG.pnts;
                origTimes = [0 (N-1)/fs];

                %% 2) FILTRADO ORIGINAL (DOS PASADAS 0.5 y 40 Hz)
                EEG=pop_eegfiltnew(EEG,'locutoff',FILT_LOW,'plotfreqz',0);
                EEG=pop_eegfiltnew(EEG,'hicutoff',FILT_HIGH,'plotfreqz',0);
                assert(EEG.pnts==N,'El filtrado cambio longitud de datos');
                filtered=EEG;
                if SAVE_FILTERED_QC
                    filtered.setname=[stem '_filtrado'];
                    pop_saveset(filtered,'filename',[stem '_filtrado.set'],'filepath',outSub);
                end

                %% 3) IDENTIFICAR CANALES (MISMOS PARAMETROS ESPACIALES ORIGINALES)
                % BurstCriterion='off' y WindowCriterion='off' evitan PRIMER ASR.
                detect = pop_clean_rawdata(filtered, ...
                    'FlatlineCriterion',FLATLINE, ...
                    'ChannelCriterion',CHAN_CORR, ...
                    'LineNoiseCriterion',LINE_NOISE, ...
                    'Highpass','off', ...
                    'BurstCriterion','off', ...
                    'WindowCriterion','off');
                goodLabels=string({detect.chanlocs.labels});
                badLabels=setdiff(originalLabels,goodLabels,'stable');
                assert(numel(goodLabels)>=MIN_GOOD_CHANNELS, ...
                    'Solo %d canales validos; revisar manualmente',numel(goodLabels));
                assert(detect.pnts==N, 'Deteccion de canales altero muestras');
                % IMPORTANTISIMO: excluir canales sobre DATOS FILTRADOS ORIGINALES
                EEG=pop_select(filtered,'channel',cellstr(goodLabels));
                assert(EEG.pnts==N,'Exclusion de canales cambio numero de muestras');
                for ci=1:numel(originalLabels)
                    CHANNELS=[CHANNELS; table(string(grp),string(subject),string(stem), ...
                        originalLabels(ci),ismember(originalLabels(ci),badLabels), ...
                        'VariableNames',{'Grupo','Sujeto','Fase','Canal','EliminadoAntesASR'})]; %#ok<AGROW>
                end

                %% 4) ASR UNICO, SOLO SOBRE CANALES VALIDOS
                EEG = pop_clean_rawdata(EEG, ...
                    'FlatlineCriterion','off','ChannelCriterion','off', ...
                    'LineNoiseCriterion','off','Highpass','off', ...
                    'BurstCriterion',BURST_CRITERION, ...
                    'WindowCriterion',WINDOW_CRITERION, ...
                    'BurstRejection','on','Distance',ASR_DISTANCE, ...
                    'WindowCriterionTolerances',WINDOW_TOLERANCES, ...
                    'MaxMem',ASR_MAXMEM);
                assert(isfield(EEG.etc,'clean_sample_mask'), ...
                    'ASR no genero clean_sample_mask');
                mask=logical(EEG.etc.clean_sample_mask(:)');
                assert(numel(mask)==N,'Mascara ASR: longitud no coincide con EDF');
                assert(nnz(mask)==EEG.pnts,'Mascara ASR no coincide con EEG.pnts');
                kept_orig_idx = find(mask); %#ok<NASGU>
                assert(all(diff(kept_orig_idx)>0),'Indices conservados no crecientes');
                if SAVE_PRE_ICA_QC
                    qc=EEG; qc.setname=[stem '_postASR_preICA'];
                    pop_saveset(qc,'filename',[stem '_postASR_preICA.set'],'filepath',outSub);
                end

                %% 5) ICA, SOLO CANALES EEG (no ECG ni canales malos)
                % No interpolar antes de ICA. Registrar rango por seguridad.
                % Rango estimado a partir de un fragmento inicial de hasta 30 s.
                % Evitar parentesis ambiguos y crear la matriz una sola vez.
                Xrank = double(EEG.data(:,1:min(EEG.pnts,round(30*fs))));
                datarank = rank(Xrank * Xrank');
                if datarank < EEG.nbchan
                    fprintf('AVISO: rango estimado %d < %d canales. ICA con PCA=%d\n', ...
                        datarank,EEG.nbchan,datarank);
                    EEG=pop_runica(EEG,'icatype','runica','extended',1,'pca',datarank, 'interrupt','off');
                else
                    EEG=pop_runica(EEG,'icatype','runica','extended',1,'interrupt','off');
                end
                EEG=pop_iclabel(EEG,'default');
                probs=EEG.etc.ic_classification.ICLabel.classifications;
                classes=string(EEG.etc.ic_classification.ICLabel.classes);
                assert(size(probs,2)>=7,'ICLabel no devolvio 7 clases');
                assert(size(probs,1)==size(EEG.icaweights,1),'Numero de IC incongruente');

                %% ECG COMO REFERENCIA EXTERNA: NO es canal ICA
                % Correlacion descriptiva, NO clasificador ni rechazo automatizado.
                % Correlacion a lag cero no capta todas las morfologias cardiacas;
                % una revision ECG-R-peaks de IC requiere un analisis posterior.
                % Se compara el ECG en los MISMOS indices originales conservados.
                ecgCorr = nan(size(probs,1),1);
                if ~isempty(ecg_orig)
                    alignedECG = ecg_orig(kept_orig_idx);
                    % Quitar componente DC, normalizar; NO modifica EEG.
                    alignedECG = alignedECG - mean(alignedECG,'omitnan');
                    icact=double(EEG.icaweights*EEG.icasphere*EEG.data(EEG.icachansind,:));
                    if std(alignedECG)>0
                        for ic=1:size(icact,1)
                            a=icact(ic,:); a=a-mean(a);
                            denom=sqrt(sum(a.^2)*sum(alignedECG.^2));
                            if denom>0
                                ecgCorr(ic)=sum(a.*alignedECG)/denom;
                            end
                        end
                    end
                end
                %% QC CARDIACO QRS: EVIDENCIA INDEPENDIENTE, NO CAMBIA RECHAZO ICA
                % Detectar R en ECG ORIGINAL, antes de compresion ASR.
                % Extraer ICs SOLO si el entorno del latido es continuo en tiempo original.
                % El ECG no se incluye en ICA y este QC NO elimina componentes.
                qrsInfo=struct('status','Sin_ECG','rOriginal',[], ...
                    'nR',0,'bpm',NaN,'nCleanR',0,'nBeatsIC',zeros(size(probs,1),1), ...
                    'lockedP2P_SD',nan(size(probs,1),1));
                if ECG_QRS_QC && ~isempty(ecg_orig)
                    try
                        qrsInfo=ecg_qrs_ic_qc(ecg_orig,fs,kept_orig_idx,icact);
                        r_peaks_original=qrsInfo.rOriginal; %#ok<NASGU>
                        save(fullfile(outSub,[stem '_ECG_QRS_QC.mat']), ...
                            'qrsInfo','r_peaks_original','-v7.3');
                    catch ME_QRS
                        qrsInfo.status=['Error_QRS: ' ME_QRS.message];
                        warning('EEG_CREA:QRS', '%s (%s)',qrsInfo.status,stem);
                    end
                end
                QRS_SUMMARY=[QRS_SUMMARY; table(string(grp),string(subject), ...
                    string(stem),string(qrsInfo.status),qrsInfo.nR,qrsInfo.bpm, ...
                    qrsInfo.nCleanR,'VariableNames',{'Grupo','Sujeto','Fase', ...
                    'EstadoQRS','PicosRDetectados','FrecuenciaCardiacaMediana_lpm', ...
                    'PicosRConVentanaEEGContinua'})]; %#ok<AGROW>
                artifactClasses = 2:6; % muscle, eye, heart, line noise, channel noise
                toRemove = any(probs(:,artifactClasses) >= ICLABEL_THRESHOLD,2);
                removedICs = find(toRemove);
                for ic=1:size(probs,1)
                    [mx,imx]=max(probs(ic,:));
                    COMPONENTS=[COMPONENTS; table(string(grp),string(subject),string(stem), ...
                        ic,classes(imx),mx,probs(ic,1),probs(ic,2),probs(ic,3), ...
                        probs(ic,4),probs(ic,5),probs(ic,6),probs(ic,7),ecgCorr(ic), ...
                        qrsInfo.nBeatsIC(ic),qrsInfo.lockedP2P_SD(ic),string(qrsInfo.status), ...
                        toRemove(ic), 'VariableNames',{'Grupo','Sujeto','Fase','IC', ...
                        'ClasePrincipal','ProbPrincipal','PBrain','PMuscle','PEye','PHeart', ...
                        'PLineNoise','PChannelNoise','POther','CorrECG_ceroLag', ...
                        'QRS_LatidosUsados','QRS_ERP_PicoPico_NormalizadoSD', ...
                        'QRS_Estado','Eliminada'})]; %#ok<AGROW>
                end
                % Evitar eliminar TODAS las ICs sin inspeccion.
                assert(nnz(~toRemove)>=2, ...
                    'ICLabel eliminaria %d/%d componentes: revisar antes de guardar', ...
                    nnz(toRemove),numel(toRemove));
                if ~isempty(removedICs), EEG=pop_subcomp(EEG,removedICs,0); end

                %% 6) INTERPOLACION ESPACIAL: al final, referencia original preservada
                % pop_interp con chanlocs originales interpola solo faltantes.
                if ~isempty(badLabels)
                    EEG=pop_interp(EEG,referenceChanlocs,'spherical');
                end
                EEG=eeg_checkset(EEG);
                assert(EEG.nbchan==numel(referenceChanlocs),'Fallo restauracion del montaje');
                assert(EEG.pnts==nnz(mask),'ICA/interpolacion alteraron numero de muestras');
                assert(~any(strcmpi(string({EEG.chanlocs.labels}),'F3')), ...
                    'F3 inesperado despues de interpolar');
                % No rereferenciar automaticamente: script historico no lo hacia.

                %% 7) TRAZABILIDAD Y FICHEROS (no interpolar HUECOS TEMPORALES)
                dm=diff([false mask false]); starts=find(dm==1); ends=find(dm==-1)-1;
                fragmentLength=ends-starts+1;
                for fr=1:numel(starts)
                    FRAGMENTS=[FRAGMENTS; table(string(grp),string(subject),string(stem), ...
                        fr,starts(fr),ends(fr),(starts(fr)-1)/fs,ends(fr)/fs, ...
                        fragmentLength(fr),fragmentLength(fr)/fs, ...
                        'VariableNames',{'Grupo','Sujeto','Fase','Fragmento', ...
                        'PrimeraMuestraOriginal','UltimaMuestraOriginal', ...
                        'InicioOriginal_s','FinOriginal_s','Muestras','Duracion_s'})]; %#ok<AGROW>
                end
                for wi=1:ceil(N/(QC_EPOCH_SECONDS*fs))
                    lo=(wi-1)*round(QC_EPOCH_SECONDS*fs)+1;
                    hi=min(wi*round(QC_EPOCH_SECONDS*fs),N);
                    local=mask(lo:hi);
                    dd=diff([false local false]); lens=find(dd==-1)-find(dd==1);
                    welch2=sum(floor(lens/(2*fs))); % ventanas de 2s sin solapamiento
                    WINDOWS=[WINDOWS; table(string(grp),string(subject),string(stem), ...
                        wi,(lo-1)/fs,hi/fs,numel(local),nnz(local), ...
                        100*nnz(local)/numel(local),numel(lens), ...
                        sum(lens(lens>=2*fs))/fs,welch2, ...
                        'VariableNames',{'Grupo','Sujeto','Fase','Ventana2min', ...
                        'InicioOriginal_s','FinOriginal_s','MuestrasOriginales', ...
                        'MuestrasConservadas','CoberturaPct','Fragmentos', ...
                        'SegundosEnFragmentosDe2s','VentanasWelch2sSinSolape'})]; %#ok<AGROW>
                end
                EEG.setname=[stem '_preprocesado'];
                EEG.etc.EEG_CREA = struct('version','v5_20261026', ...
                    'sourceEDF',srcEDF,'group',grp,'subject',subject,'phase',stem, ...
                    'originalSampleCount',N,'sampleRate',fs, ...
                    'badChannels',{cellstr(badLabels)}, ...
                    'originalChanlocs',{referenceChanlocs}, ...
                    'originalSampleIndexAfterASR',kept_orig_idx, ...
                    'parameters',struct('filterHz',[FILT_LOW FILT_HIGH], ...
                    'FlatlineCriterion',FLATLINE,'ChannelCriterion',CHAN_CORR, ...
                    'LineNoiseCriterion',LINE_NOISE,'BurstCriterion',BURST_CRITERION, ...
                    'WindowCriterion',WINDOW_CRITERION, ...
                    'WindowCriterionTolerances',WINDOW_TOLERANCES, ...
                    'BurstRejection','on','ICLabelThreshold',ICLABEL_THRESHOLD));
                % Guardar mapa completo, no inferir tiempo desde muestras comprimidas.
                maskOrig = mask; %#ok<NASGU>
                original_sample_index = kept_orig_idx; %#ok<NASGU>
                original_srate=fs; %#ok<NASGU>
                original_n_samples=N; %#ok<NASGU>
                original_time_span_s=origTimes; %#ok<NASGU>
                save(fullfile(outSub,[stem '_MAPA_TEMPORAL.mat']), ...
                    'maskOrig','original_sample_index','original_srate', ...
                    'original_n_samples','original_time_span_s','-v7.3');
                if SAVE_ECG_REFERENCE && ~isempty(ecg_orig)
                    ecg_original=ecg_orig; %#ok<NASGU>
                    save(fullfile(outSub,[stem '_ECG_REFERENCIA.mat']), ...
                        'ecg_original','original_srate','original_n_samples','-v7.3');
                end
                EEG=pop_saveset(EEG,'filename',[stem '_preprocesado.set'], ...
                    'filepath',outSub);
                if SAVE_FULL_MAT
                    save(fullfile(outSub,['EEG_preprocess_' stem '.mat']),'EEG','-v7.3');
                end
                %% 8) VERIFICAR QUE LO GUARDADO PUEDE RECARGARSE SIN PERDER TRAZABILIDAD
                if VERIFY_SAVED_FILES
                    verify_saved_crea(outSub,stem,N,fs,mask,kept_orig_idx, ...
                        referenceChanlocs,EEG,SAVE_FULL_MAT,REQUIRE_NO_F3);
                end
                %% 9) PSD QC: comparar filtrado original vs EEG final
                % Se grafican todos los 19 canales, mismos colores y eje Y.
                % DESPUES: Welch por tramo temporal continuo; JAMAS cruzar boundaries.
                % No es un analisis confirmatorio ni sustituye PSD/FOOOF definitivo.
                if SAVE_PSD_QC
                    try
                        qcinfo = plot_psd_qc_crea(filtered, EEG, mask, referenceChanlocs, ...
                            fs, stem, outSub, badLabels, PSD_QC_WINDOW_SEC, PSD_QC_YLIM_COMMON);
                        fprintf('PSD_QC %s: %d ventanas antes, %d despues (%s)\n', ...
                            stem,qcinfo.nBefore,qcinfo.nAfter,qcinfo.file);
                    catch ME_QC
                        warning('EEG_CREA:PSDQCFallo', ...
                            'Preprocesamiento guardado, pero fallo grafico PSD %s: %s',stem,ME_QC.message);
                        INCIDENTS=[INCIDENTS; table(string(grp),string(subject),string(stem), ...
                            string(['PSDQCFallo: ' ME_QC.message]),'VariableNames', ...
                            {'Grupo','Sujeto','Fase','Incidencia'})]; %#ok<AGROW>
                    end
                end
                % Marcar como completo solo cuando han terminado guardado,
                % verificacion y QC (si falla PSD QC queda incidencia, no completo).
                qcPNG = fullfile(outSub,[stem '_PSD_QC_ANTES_DESPUES.png']);
                if VERIFY_SAVED_FILES && (~SAVE_PSD_QC || isfile(qcPNG))
                    fid=fopen(completeMarker,'w');
                    assert(fid>0,'No se pudo escribir marcador de integridad');
                    fprintf(fid,'VERIFICADO %s\n',datestr(now,31)); fclose(fid);
                end
                SUMMARY=[SUMMARY; table(string(grp),string(subject),string(stem), ...
                    N,nnz(mask),N-nnz(mask),100*(1-nnz(mask)/N), ...
                    numel(originalLabels),numel(goodLabels),numel(badLabels), ...
                    string(strjoin(badLabels,',')),size(probs,1),nnz(toRemove), ...
                    string(mat2str(removedICs(:)')),sum(fragmentLength>=2*fs), ...
                    max([0 fragmentLength])/fs,logical(~isempty(ecg_orig)), ...
                    'VariableNames',{'Grupo','Sujeto','Fase','MuestrasOriginales', ...
                    'MuestrasConservadas','MuestrasRechazadas','RechazoPct', ...
                    'CanalesOriginales','CanalesParaASR_ICA','NumCanalesMalos', ...
                    'CanalesMalos','NumIC','NumICEliminadas','ICEliminadas', ...
                    'FragmentosAlMenos2s','TramoContinuoMax_s','ECGDisponible'})]; %#ok<AGROW>
                fprintf('OK %s: rechazadas %.2f%%; canales malos {%s}; IC eliminadas %d/%d\n', ...
                    stem,100*(1-nnz(mask)/N),strjoin(badLabels,','), ...
                    nnz(toRemove),size(probs,1));
            catch ME
                fprintf(2,'ERROR en %s: %s\n',stem,ME.getReport('basic'));
                INCIDENTS=[INCIDENTS; table(string(grp),string(subject),string(stem), ...
                    string(ME.message),'VariableNames',{'Grupo','Sujeto','Fase','Incidencia'})]; %#ok<AGROW>
            end
            % Guardados incrementales: si MATLAB se interrumpe, ya hay QC.
            write_audit_files(QCDIR,SUMMARY,CHANNELS,COMPONENTS,QRS_SUMMARY,WINDOWS,FRAGMENTS,INCIDENTS,MIN_WELCH_WINDOWS_2MIN,PSD_QC_WINDOW_SEC);
        end
    end
end
write_audit_files(QCDIR,SUMMARY,CHANNELS,COMPONENTS,QRS_SUMMARY,WINDOWS,FRAGMENTS,INCIDENTS,MIN_WELCH_WINDOWS_2MIN,PSD_QC_WINDOW_SEC);
fprintf('\nTERMINADO. Informe: %s\n',QCDIR);


%% VALIDACION POST-GUARDADO: INTEGRIDAD DE MUESTRAS, ETIQUETAS Y HUECOS
function verify_saved_crea(outSub,stem,N,fs,mask,idx,chanlocs,EEG_expected,fullMat,requireNoF3)
    expectedMask=logical(mask(:)');
    expectedIdx=double(idx(:)');
    assert(numel(expectedMask)==N && nnz(expectedMask)==numel(expectedIdx), ...
        'Mascara original inconsistente antes de guardar');
    assert(isequal(expectedIdx,double(find(expectedMask))), ...
        'Indices originales no corresponden a la mascara');
    assert(numel(expectedIdx)==EEG_expected.pnts, ...
        'Conteo de muestras incorrecto antes del guardado');
    % Validacion de todos los fragmentos: dentro de ellos no puede haber salto.
    cuts=[1 find(diff(expectedIdx)>1)+1 numel(expectedIdx)+1];
    for q=1:numel(cuts)-1
        v=expectedIdx(cuts(q):cuts(q+1)-1);
        assert(all(diff(v)==1),'Fragmento atraviesa discontinuidad temporal');
    end
    mapfile=fullfile(outSub,[stem '_MAPA_TEMPORAL.mat']);
    M=load(mapfile,'maskOrig','original_sample_index', ...
        'original_n_samples','original_srate');
    assert(M.original_n_samples==N && abs(M.original_srate-fs)<1e-10, ...
        'Metadatos temporales incorrectos al recargar el mapa');
    assert(isequal(logical(M.maskOrig(:)'),expectedMask), ...
        'Mascara distinta al recargar el mapa');
    assert(isequal(double(M.original_sample_index(:)'),expectedIdx), ...
        'Indices originales distintos al recargar el mapa');
    setfile=fullfile(outSub,[stem '_preprocesado.set']);
    R=pop_loadset('filename',[stem '_preprocesado.set'],'filepath',outSub);
    assert(R.pnts==numel(expectedIdx) && R.nbchan==numel(chanlocs), ...
        'Dimension incorrecta al recargar .set');
    assert(abs(R.srate-fs)<1e-10,'Fs incorrecta al recargar .set');
    assert(isfield(R.etc,'clean_sample_mask') && ...
        isequal(logical(R.etc.clean_sample_mask(:)'),expectedMask), ...
        'clean_sample_mask perdida o alterada en .set');
    assert(isfield(R.etc,'EEG_CREA') && ...
        isfield(R.etc.EEG_CREA,'originalSampleIndexAfterASR') && ...
        isequal(double(R.etc.EEG_CREA.originalSampleIndexAfterASR(:)'),expectedIdx), ...
        'Correspondencia temporal alterada en .set');
    assert(isequal(string({R.chanlocs.labels}),string({chanlocs.labels})), ...
        'Montaje/canales alterados al recargar .set');
    if requireNoF3
        assert(~any(strcmpi(string({R.chanlocs.labels}),'F3')), ...
            'F3 inesperado en el .set recargado');
    end
    % Igualdad de datos tras escritura/lectura (tolerancia para precision single).
    x=double(R.data); y=double(EEG_expected.data);
    assert(isequal(size(x),size(y)),'Forma de datos alterada al recargar .set');
    deviation=max(abs(x(:)-y(:)));
    assert(isfinite(deviation) && deviation<=1e-4, ...
        'Datos del .set alterados (max |diferencia|=%.9g)',deviation);
    if fullMat
        D=load(fullfile(outSub,['EEG_preprocess_' stem '.mat']),'EEG');
        Q=D.EEG;
        assert(Q.pnts==R.pnts && Q.nbchan==R.nbchan, ...
            'Dimensiones .mat distintas de .set');
        assert(isfield(Q.etc,'clean_sample_mask') && ...
            isequal(logical(Q.etc.clean_sample_mask(:)'),expectedMask), ...
            'Mascara alterada al recargar .mat');
        assert(isfield(Q.etc,'EEG_CREA') && ...
            isequal(double(Q.etc.EEG_CREA.originalSampleIndexAfterASR(:)'),expectedIdx), ...
            'Indices alterados al recargar .mat');
        assert(isequal(string({Q.chanlocs.labels}),string({R.chanlocs.labels})), ...
            'Canales distintos entre .mat y .set');
    end
    % Comprobar marcadores de discontinuidad cuando existan saltos originales.
    gaps=sum(diff(expectedIdx)>1);
    if gaps>0
        types=string({R.event.type});
        nb=sum(strcmpi(types,'boundary'));
        assert(nb>=gaps, ...
            'Al recargar .set: %d saltos originales, pero solo %d boundaries',gaps,nb);
    end
    fprintf('QC RECARGA OK %s: %d muestras, %d saltos, integridad .set/.mat/mapa verificada\n', ...
        stem,numel(expectedIdx),gaps);
end

%% FUNCIONES LOCALES
function write_audit_files(folder,summary,channels,components,qrsSummary,windows,fragments,incidents,minWelch,welchWindowSec)
    % WELCH_2MIN es una AUDITORIA DE DISPONIBILIDAD, no PSD estimada.
    % Cada ventana de Welch debe estar completamente dentro de un fragmento
    % conservado, y cada fragmento dentro del intervalo de 2 minutos original.
    % El conteo 'VentanasWelch2sSinSolape' se calculo separadamente para cada
    % intervalo; no se permite cruzar un boundary ni un limite de intervalo.
    welch = windows;
    if ~isempty(welch)
        welch.VentanaWelch_s = repmat(welchWindowSec,height(welch),1);
        welch.VentanasMinimasProvisionales = repmat(minWelch,height(welch),1);
        welch.SegundosEfectivosWelch = welch.VentanasWelch2sSinSolape * welchWindowSec;
        welch.PSDCalculableProvisional = welch.VentanasWelch2sSinSolape >= minWelch;
        welch.EstadoPSD = repmat("NaN_insuficientes_ventanas",height(welch),1);
        welch.EstadoPSD(welch.PSDCalculableProvisional) = "Calcular_PSD";
        % El umbral no excluye sujetos ni descarta datos preprocesados.
        % PSD/FOOOF definitivos aplicaran esta regla y escribiran NaN donde toque.
    end
    tabs={summary,channels,components,qrsSummary,windows,welch,fragments,incidents};
    names={'RESUMEN','CANALES','COMPONENTES','ECG_QRS','COBERTURA_2MIN','WELCH_2MIN','FRAGMENTOS','INCIDENCIAS'};
    for z=1:numel(tabs)
        if ~isempty(tabs{z})
            writetable(tabs{z},fullfile(folder,[names{z} '.csv']));
        end
    end
    % Excel es una copia cómoda de CSV, no sustituye MAT de trazabilidad.
    xls=fullfile(folder,'AUDITORIA_PREPROCESAMIENTO.xlsx');
    if isfile(xls), delete(xls); end
    for z=1:numel(tabs)
        if ~isempty(tabs{z})
            writetable(tabs{z},xls,'Sheet',names{z});
        end
    end
end


%% PSD QC: NO cruzar discontinuidades del ASR
function info = plot_psd_qc_crea(filtered, cleaned, mask, allChanlocs, fs, stem, dest, badLabels, epochSeconds, commonY)
    labels = string({allChanlocs.labels});
    % La interpolacion debe devolver las etiquetas y el numero originales.
    cleanLabels = string({cleaned.chanlocs.labels});
    [exists, reorder] = ismember(labels,cleanLabels);
    assert(all(exists),'Faltan canales al graficar PSD QC');
    Xbefore=double(filtered.data);
    origLabels=string({filtered.chanlocs.labels});
    [allPresent, orderBefore] = ismember(labels,origLabels);
    assert(all(allPresent),'Faltan canales en EEG filtrado');
    Xbefore=Xbefore(orderBefore,:);
    Xafter=double(cleaned.data(reorder,:));
    nSamp=round(epochSeconds*fs);
    assert(nSamp>=32,'Ventana de Welch demasiado corta');
    % Del mismo tamaño de ventana y nfft en ambas condiciones.
    win=hamming(nSamp);
    nfft=nSamp;
    [pBefore,f,nBefore]=psd_from_intact_fragments(Xbefore, ...
        [1 size(Xbefore,2)],win,nfft,fs);
    idx=find(logical(mask(:)'));
    assert(numel(idx)==size(Xafter,2),'Mascara inconsistente para PSD');
    if isempty(idx), error('Sin muestras despues de ASR'); end
    % Rango de indices en EEG FINAL por cada tramo ORIGINAL contiguo.
    starts=[1 find(diff(idx)>1)+1];
    ends=[starts(2:end)-1 numel(idx)];
    [pAfter,fa,nAfter]=psd_from_intact_fragments(Xafter, ...
        [starts(:) ends(:)],win,nfft,fs);
    assert(isequal(f,fa),'Frecuencias no comparables');
    assert(nBefore>0 && nAfter>0,'Sin ventanas Welch completas de 2 s');
    dBbefore=10*log10(max(pBefore,realmin));
    dBafter=10*log10(max(pAfter,realmin));
    keep=f>=0.5 & f<=40;
    colors=lines(numel(labels));
    h=figure('Visible','off','Color','w','Position',[60 80 1600 690]);
    onCleanupFigure=onCleanup(@()close(h)); %#ok<NASGU>
    ax1=subplot(1,2,1,'Parent',h); hold(ax1,'on');
    ax2=subplot(1,2,2,'Parent',h); hold(ax2,'on');
    for c=1:numel(labels)
        plot(ax1,f(keep),dBbefore(c,keep),'LineWidth',0.9,'Color',colors(c,:));
        plot(ax2,f(keep),dBafter(c,keep),'LineWidth',0.9,'Color',colors(c,:));
    end
    title(ax1,sprintf('ANTES: filtrado 0.5-40 Hz | %d ventanas',nBefore),'Interpreter','none');
    title(ax2,sprintf('DESPUES: ASR + ICA + interpolacion | %d ventanas',nAfter),'Interpreter','none');
    for ax=[ax1 ax2]
        xlabel(ax,'Frecuencia (Hz)');
        ylabel(ax,'PSD (10 log_{10}[\muV^2/Hz])');
        xlim(ax,[0.5 40]); grid(ax,'on'); box(ax,'on');
    end
    if commonY
        v=[dBbefore(:,keep) dBafter(:,keep)];
        v=v(isfinite(v));
        if ~isempty(v)
            yl=[floor(min(v(:))) ceil(max(v(:)))];
            if yl(1)==yl(2), yl=yl+[-1 1]; end
            ylim(ax1,yl); ylim(ax2,yl);
        end
    end
    lg=legend(ax2,cellstr(labels),'Location','eastoutside');
    lg.FontSize=8;
    titleStr=sprintf('%s | PSD QC pre/post | canales interpolados: %s', ...
        stem,strjoin(badLabels,', '));
    annotation(h,'textbox',[0.02 0.94 0.96 0.05], 'String',titleStr, ...
        'LineStyle','none','FontSize',14,'FontWeight','bold', ...
        'HorizontalAlignment','center','Interpreter','none');
    annotation(h,'textbox',[0.02 0.005 0.96 0.055], ...
        'String','QC descriptivo: Welch 2 s Hamming, sin solape; DESPUES calculado dentro de tramos originales continuos (sin cruzar boundaries). La interpolacion modifica los canales marcados.', ...
        'LineStyle','none','FontSize',9,'HorizontalAlignment','center');
    file=fullfile(dest,[stem '_PSD_QC_ANTES_DESPUES.png']);
    print(h,file,'-dpng','-r160');
    info=struct('file',file,'nBefore',nBefore,'nAfter',nAfter);
    % Guardar valores reproducibles + mapeo de canales para revisarlos.
    save(fullfile(dest,[stem '_PSD_QC_VALORES.mat']), ...
        'pBefore','pAfter','f','labels','badLabels','nBefore','nAfter','-v7.3');
end

function [psd,f,nWin]=psd_from_intact_fragments(X,ranges,window,nfft,fs)
    % Promedio de TODAS las ventanas de longitud fija; ninguna cruza huecos.
    L=numel(window);
    nCh=size(X,1);
    f=(0:nfft/2)'*fs/nfft;
    psd=zeros(nCh,numel(f));
    nWin=0;
    for r=1:size(ranges,1)
        a=ranges(r,1); b=ranges(r,2);
        n=floor((b-a+1)/L);
        if n<1, continue; end
        % pwelch sobre cada fragmento CONTINUO; pesos segun ventanas completas.
        seg=X(:,a:a+n*L-1);
        for c=1:nCh
            [p,fr]=pwelch(seg(c,:),window,0,nfft,fs);
            assert(numel(p)==numel(f),'Eje de frecuencia inesperado');
            if c==1, f=fr; end
            psd(c,:)=psd(c,:)+n*p(:)';
        end
        nWin=nWin+n;
    end
    if nWin>0, psd=psd/nWin; else, psd(:)=NaN; end
end


%% ECG-QRS: MEDIDA DESCRIPTIVA DE ACTIVIDAD DE IC SINCRONIZADA CON LATIDOS
function out=ecg_qrs_ic_qc(ecg_orig,fs,idx,icact)
    % Requiere Signal Processing Toolbox (butter, filtfilt, findpeaks).
    % NO se aplican filtros al EEG; ECG se analiza en su propio eje original.
    nIC=size(icact,1); N=numel(ecg_orig);
    out=struct('status','', 'rOriginal',[], 'nR',0, 'bpm',NaN, ...
        'nCleanR',0,'nBeatsIC',zeros(nIC,1), ...
        'lockedP2P_SD',nan(nIC,1));
    ecg=double(ecg_orig(:)');
    if N<fs*30 || any(~isfinite(ecg)) || std(ecg)==0
        out.status='ECG_invalido_corto_o_no_finito'; return;
    end
    [b,a]=butter(3,[5 25]/(fs/2),'bandpass');
    xf=filtfilt(b,a,ecg);
    % ECG bipolar: buscar prominencias en ambas polaridades; elegir la que
    % produce mas picos de frecuencia cardiaca plausible.
    minDist=round(.4*fs);
    noise=median(abs(xf-median(xf)))*1.4826;
    prominence=max(3*noise,eps);
    [~,pPlus]=findpeaks(xf,'MinPeakDistance',minDist, ...
        'MinPeakProminence',prominence);
    [~,pMinus]=findpeaks(-xf,'MinPeakDistance',minDist, ...
        'MinPeakProminence',prominence);
    opts={pPlus,pMinus}; quality=zeros(1,2);
    for q=1:2
        rr=diff(opts{q})/fs;
        if isempty(rr), continue; end
        bpm=60/median(rr);
        quality(q)=numel(opts{q}) * double(bpm>=40 && bpm<=180);
    end
    if all(quality==0)
        out.status='Picos_R_no_fiables_revisar_ECG'; return;
    end
    [~,best]=max(quality); r=opts{best};
    rr=diff(r)/fs;
    % Evitar aceptar ritmos muy irregulares o deteccion evidentemente fallida.
    medianBPM=60/median(rr);
    out.rOriginal=r(:)'; out.nR=numel(r); out.bpm=medianBPM;
    out.status='QRS_detectado_QC_pendiente_revision_visual';
    % Indexacion inversa (0 para muestras eliminadas). No cruzar huecos.
    idx=double(idx(:)');
    map=zeros(1,N,'uint32'); map(idx)=uint32(1:numel(idx));
    pre=round(.2*fs); post=round(.3*fs);
    lags=-pre:post; windows=zeros(nIC,numel(lags)); count=0;
    for j=1:numel(r)
        orig=r(j)+lags;
        if orig(1)<1 || orig(end)>N, continue; end
        loc=double(map(orig));
        if any(loc==0) || any(diff(loc)~=1), continue; end
        windows=windows+double(icact(:,loc));
        count=count+1;
    end
    out.nCleanR=count; out.nBeatsIC(:)=count;
    if count<10
        out.status='Menos_de_10_latidos_R_con_entorno_EEG_continuo';return;
    end
    erp=windows/count;
    % Quitar linea base antes de R, sin modificar EEG ni las IC originales.
    baseline=mean(erp(:,1:max(1,round(.1*fs))),2);
    erp=erp-baseline;
    sd=std(double(icact),0,2);
    sd(sd==0)=NaN;
    out.lockedP2P_SD=(max(erp,[],2)-min(erp,[],2))./sd;
    % Este indice es DESCRIPTIVO: puede reflejar actividad neural
    % cardiaca-evocada; no diagnostica artefacto y no elimina IC.
end
