%% EEG_CREA_FOOOF_V9_ORIGINAL_GRUPO_INCLUIDO.m
% EEG_CREA | FOOOF 1.x (implementacion Python original llamada desde MATLAB)
% Entrada PSD: PSD_LONG_2MIN.csv ORIGINAL sin limpieza IQR (PSD lineal).
% La columna Grupo (1 o 2) ya figura en EL MISMO CSV original.
% No se usa ningun archivo auxiliar ni PSD procesada por IQR.
% 41 sujetos x 19 electrodos x 21 momentos x 80 frecuencias = 1.308.720 filas.
% Se ajusta SOLO 2:0.5:40 Hz (77 frecuencias); NO interpolar NaN.
% IMPORTANTE: NO renumerar ni compactar los momentos al excluir espectros.
% M1 PRE_OA, M2:M16 tSMS_1:tSMS_15, M17:M21 POST_1:POST_5.
% Requisitos: MATLAB R2019b+ (pyenv), Python compatible, pip install fooof==1.1.1
%           y numpy. Se requiere la API fooof.FOOOF (NO specparam 2.x).
% Salidas: 01_RESULTADOS_FOOOF.xlsx, 02_TODOS_LOS_PICOS_FOOOF.xlsx, 03_AUDITORIA_FOOOF.xlsx
%          figuras de QC en FIGURAS_FOOOF. No genera CSV.
% El modelo se ajusta sobre potencia LINEAL: FOOOF realiza el log internamente.
% R2 principal = fm.r_squared_ de FOOOF; R2_Corr se guarda como QC auxiliar
% (correlacion cuadrada entre espectro log10 y modelo completo log10).
% No se aplica umbral de R2 para excluir parametros automaticamente.
% V9: Grupo se lee del CSV original; sin segundo CSV. R2 oficial de FOOOF.
% Delta FOOOF refleja solo picos con centro >=2 Hz (delta parcial).
% Nota: la validacion MATLAB-Python definitiva requiere ejecutar en tu equipo.

clearvars; close all; clc;
%% ================= CONFIGURACION =================
BASE = 'C:\Users\evadi\OneDrive\Escritorio\FENNSI\CREATIVIDAD_DEFINITIVOS_EEG\Database_Resultados\PSD';
ENTRADA = fullfile(BASE,'PSD_LONG_2MIN.csv');
SALIDA = 'C:\Users\evadi\OneDrive\Escritorio\FENNSI\CREATIVIDAD_DEFINITIVOS_EEG\Database_Resultados\APERIODICO_ORIGINAL_V9';
% Ajusta ENTRADA si la ubicacion del CSV original con Grupo es distinta.
RANGO = [2 40];
FREQ_80 = (0.5:0.5:40)';
FREQ_FIT = (2:0.5:40)';
ANCHO_PICOS = [1 12];  % 1 Hz = 2 bins para resolucion de 0.5 Hz
MAX_PICOS = 6;      % limite elegido; el default de FOOOF es ilimitado
MIN_ALTURA = 0;       % default
UMBRAL_PICOS = 2;    % default
MODO_AP = 'fixed';
HACER_FIGURAS = true;
N_EJEMPLOS = 6;       % ejemplos de ajustes individuales, no sinteticos

if ~isfile(ENTRADA)
    [archivoSeleccionado, carpetaSeleccionada] = uigetfile({'*.csv','Archivos CSV (*.csv)'}, ...
        'Selecciona PSD_LONG_2MIN.csv ORIGINAL (sin IQR)');
    assert(~isequal(archivoSeleccionado,0), ...
        'No se selecciono el CSV de PSD originales.');
    ENTRADA = fullfile(carpetaSeleccionada,archivoSeleccionado);
end
assert(isfile(ENTRADA),'No existe archivo entrada: %s',ENTRADA);
if ~exist(SALIDA,'dir'), mkdir(SALIDA); end
CARPETA_FIG = fullfile(SALIDA,'FIGURAS_FOOOF');
if HACER_FIGURAS && ~exist(CARPETA_FIG,'dir'),mkdir(CARPETA_FIG);end
OUT_RESULTADOS = fullfile(SALIDA,'01_RESULTADOS_FOOOF.xlsx');
OUT_PICOS = fullfile(SALIDA,'02_TODOS_LOS_PICOS_FOOOF.xlsx');
OUT_AUDITORIA = fullfile(SALIDA,'03_AUDITORIA_FOOOF.xlsx');
% Evitar ejecutar miles de ajustes si existen los Excel de una corrida anterior.
assert(~isfile(OUT_RESULTADOS) && ~isfile(OUT_PICOS) && ~isfile(OUT_AUDITORIA), ...
    'Ya existen archivos FOOOF en SALIDA. Cambia la carpeta de salida antes de ejecutar.');

%% ================= VERSIONES / MOTOR =================
pe = pyenv;
fprintf('Python: %s\n',char(pe.Version));
try
    np = py.importlib.import_module('numpy');
    fooofpy = py.importlib.import_module('fooof');
    versionFooof = char(py.str(py.getattr(fooofpy,'__version__')));
    fprintf('FOOOF Python: %s\n',versionFooof);
    assert(startsWith(versionFooof,'1.'), ...
        'Este script requiere fooof 1.x. Version instalada: %s',versionFooof);
catch ME
    error('No se pudo conectar con FOOOF 1.x desde MATLAB: %s',ME.message);
end

%% ================= LECTURA / VALIDACION =================
fprintf('Leyendo PSD ORIGINALES sin limpieza IQR...\n');
T = readtable(ENTRADA,'TextType','string');
req = {'Sujeto','Grupo','Electrodo','Momento','Fase','Frecuencia_Hz','PSD_lineal'};
assert(all(ismember(req,T.Properties.VariableNames)),'Faltan columnas requeridas en PSD original.');
assert(height(T)==41*19*21*80,'Numero inesperado de filas: %d',height(T));
T.Sujeto = string(T.Sujeto);
% Comprobar que el CSV es el ORIGINAL (sin IQR), aunque contenga Grupo.
% Esta prueba se basa en los recuentos exactos de la auditoria original.
assert(nnz(isnan(T.PSD_lineal))==24320 && ...
    nnz(isfinite(T.PSD_lineal) & T.PSD_lineal>0)==1284400, ...
    ['La entrada no reproduce la PSD original auditada: se esperan ' ...
    '24.320 NaN y 1.284.400 valores positivos validos. No usar CSV IQR.']);
T.Grupo = str2double(string(T.Grupo));
assert(all(isfinite(T.Grupo) & ismember(T.Grupo,[1 2])), ...
    'Grupo debe contener solo 1 o 2 y no tener ausentes.');
[gs,~,gi] = unique(T.Sujeto);
for ii=1:numel(gs)
    assert(numel(unique(T.Grupo(gi==ii)))==1, ...
        'Grupo contradictorio en el CSV para el sujeto %s',char(gs(ii)));
end
fprintf('Grupos del CSV original: %d sujetos G1 y %d sujetos G2.\n', ...
    nnz(accumarray(gi,T.Grupo,[],@min)==1), ...
    nnz(accumarray(gi,T.Grupo,[],@min)==2));
assert(nnz(accumarray(gi,T.Grupo,[],@min)==1)==20 && ...
    nnz(accumarray(gi,T.Grupo,[],@min)==2)==21, ...
    'Numero de sujetos por grupo distinto del auditado (20 y 21).');
assert(numel(unique(T.Sujeto))==41,'No hay exactamente 41 sujetos.');
assert(all(~ismissing(T.Sujeto)),'Hay sujetos sin identificador.');
assert(all(ismember(T.Electrodo,1:19)) && numel(unique(T.Electrodo))==19, 'Electrodos incorrectos');
assert(all(ismember(T.Momento,1:21)) && numel(unique(T.Momento))==21, 'Momentos incorrectos');
assert(all(ismember(T.Frecuencia_Hz,FREQ_80)),'Frecuencias fuera de la malla esperada.');
assert(all(isfinite(T.Momento)) && all(isfinite(T.Electrodo)), ...
    'Momento o electrodo no finito.');
assert(all(isnan(T.PSD_lineal) | (isfinite(T.PSD_lineal)&T.PSD_lineal>0)), ...
    'Existen PSD no finitas distintas de NaN o no positivas.');

% Una sola malla, orden explicito, sin deducir la cronologia de las filas.
T = sortrows(T,{'Sujeto','Electrodo','Momento','Frecuencia_Hz'});
subjects=unique(T.Sujeto,'sorted');
N=41*19*21;
assert(all(reshape(T.Frecuencia_Hz,80,[])==repmat(FREQ_80,1,N),'all'), ...
    'Frecuencias duplicadas/ausentes o malla mal ordenada en algun espectro.');
PSD80 = reshape(T.PSD_lineal,80,N);
ID_S = reshape(T.Sujeto,80,N);ID_E = reshape(T.Electrodo,80,N);ID_M=reshape(T.Momento,80,N);
ID_G=reshape(T.Grupo,80,N);
assert(all(ID_S==ID_S(1,:),'all') && all(ID_E==ID_E(1,:),'all') && ...
    all(ID_M==ID_M(1,:),'all') && all(ID_G==ID_G(1,:),'all'), ...
    'La agrupacion sujeto-electrodo-momento-grupo no es homogenea.');
Sujeto=ID_S(1,:)';Electrodo=ID_E(1,:)';Momento=ID_M(1,:)';
Grupo=ID_G(1,:)';
assert(numel(unique(strcat(Sujeto,"_",string(Electrodo),"_",string(Momento))))==N, ...
    'Clave sujeto-electrodo-momento repetida o ausente.');
for s=1:numel(subjects)
    sel=Sujeto==subjects(s);
    assert(nnz(sel)==19*21,'Espectros incompletos en sujeto %s',subjects(s));
end

canales={'Fp1','Fp2','F4','F7','F8','Fz','C3','C4','Cz','T3','T4','T5','T6','P3','P4','Pz','O1','O2','Oz'};
ElectrodoNombre=string(canales(Electrodo));
ElectrodoNombre=ElectrodoNombre(:); % obligatorio: 16359 x 1, no 1 x 16359
Fase=repmat("",N,1);NombreMomento=repmat("",N,1);
Fase(Momento==1)="PRE_OA"; NombreMomento(Momento==1)="PRE_OA";
iT=Momento>=2 & Momento<=16;Fase(iT)="tSMS";NombreMomento(iT)="tSMS_"+string(Momento(iT)-1);
iP=Momento>=17 & Momento<=21;Fase(iP)="POST";NombreMomento(iP)="POST_"+string(Momento(iP)-16);
assert(all(reshape(string(T.Fase),80,N)==repmat(Fase',80,1),'all'), ...
    'Inconsistencia entre codigo Momento y Fase de origen.');

% Grupo validado directamente en el CSV ORIGINAL; no sobrescribirlo.
PSDfit=PSD80(FREQ_80>=RANGO(1)&FREQ_80<=RANGO(2),:);
assert(size(PSDfit,1)==77 && all(FREQ_FIT==FREQ_80(FREQ_80>=2 & FREQ_80<=40)));
NValidas=sum(isfinite(PSDfit)&PSDfit>0,1)';
NNaN_Rango=77-NValidas;
NNaN_Global=sum(~isfinite(PSD80),1)';
EsCompleto=(NValidas==77);

%% ================= INICIALIZAR SALIDAS =================
Offset=nan(N,1);Exponent=nan(N,1);R2=nan(N,1);R2_Corr=nan(N,1);
ErrorFOOOF=nan(N,1);NumeroPicos=nan(N,1);
Estado=repmat("SIN_AJUSTAR",N,1);MensajeError=repmat("",N,1);
% Almacenar las curvas de ejemplos solo (sin 16k x 77 en Excel).
Ejemplos=struct('Sujeto',{},'Electrodo',{},'Momento',{},'Frec',{},'Original',{},'Full',{},'AP',{});
Picos = cell(N,1);

%% ================= EJECUCION FOOOF =================
fprintf('Espectros totales: %d; completos 2-40 Hz: %d; incompletos: %d\n', ...
    N,nnz(EsCompleto),nnz(~EsCompleto));
% Prueba real ANTES del procesamiento masivo: ajusta un espectro
% y verifica el acceso a todos los atributos Python necesarios.
idxPrueba=find(EsCompleto,1);
assert(~isempty(idxPrueba),'No existe ningun espectro completo para ajustar.');
fprintf('Probando FOOOF con el espectro %d antes de procesar el lote...\n',idxPrueba);
try
    fmTest=fooofpy.FOOOF(py.list(ANCHO_PICOS),MAX_PICOS,MIN_ALTURA, ...
        UMBRAL_PICOS,MODO_AP,false);
    fmTest.fit(np.array(FREQ_FIT'),np.array(PSDfit(:,idxPrueba)'),py.list(RANGO));
    apTest=double(py.array.array('d',fmTest.aperiodic_params_));
    curveTest=double(py.array.array('d',fmTest.fooofed_spectrum_));
    ppTest=np.asarray(fmTest.peak_params_);
    if int64(py.len(ppTest)) == 0
        pkTest=zeros(0,3);
    else
        vTest=double(py.array.array('d',py.numpy.nditer(ppTest)));
        assert(mod(numel(vTest),3)==0,'Prueba inicial: formato de picos inesperado.');
        pkTest=reshape(vTest,3,[])';
    end
    aperCurveTest=double(py.array.array('d',py.getattr(fmTest,'_ap_fit')));
    assert(numel(apTest)==2 && all(isfinite(apTest)) && ...
        numel(curveTest)==77 && all(isfinite(curveTest)) && ...
        mod(numel(pkTest),3)==0 && numel(aperCurveTest)==77, ...
        'El motor Python ha devuelto parametros o curvas incompatibles.');
    clear fmTest apTest curveTest pkTest aperCurveTest ppTest
    fprintf('Prueba FOOOF OK. Iniciando procesamiento completo.\n');
catch ME
    error('PRUEBA INICIAL FOOOF FALLIDA. No se lanza el lote: %s',ME.message);
end

for k=1:N
    if ~EsCompleto(k)
        if NValidas(k)==0, Estado(k)="SIN_PSD_EN_RANGO";
        else, Estado(k)="ESPECTRO_INCOMPLETO"; end
        continue
    end
    try
        % Objeto nuevo por espectro: sin dependencia entre ajustes consecutivos.
        fm=fooofpy.FOOOF(py.list(ANCHO_PICOS),MAX_PICOS,MIN_ALTURA, ...
            UMBRAL_PICOS,MODO_AP,false);
        fm.fit(np.array(FREQ_FIT'),np.array(PSDfit(:,k)'),py.list(RANGO));
        ap=double(py.array.array('d',fm.aperiodic_params_));
        assert(numel(ap)==2 && all(isfinite(ap)), 'Parametros aperiodicos no validos');
        Offset(k)=ap(1);Exponent(k)=ap(2);
        ErrorFOOOF(k)=double(fm.error_);
        yy=double(py.array.array('d',fm.power_spectrum));
        yh=double(py.array.array('d',fm.fooofed_spectrum_));
        assert(numel(yy)==77 && numel(yh)==77 && all(isfinite(yy)) && all(isfinite(yh)), ...
            'La curva ajustada no contiene los 77 puntos finitos.');
        % R2_Corr es SOLO una metrica auxiliar: si es indefinida,
        % no debe invalidar un ajuste aperiodico correcto.
        if std(yy(:))>0 && std(yh(:))>0
            cc=corrcoef(yy(:),yh(:));
            if isequal(size(cc),[2 2]) && all(isfinite(cc),'all')
                R2_Corr(k)=cc(1,2)^2;
            end
        end
        R2(k)=double(fm.r_squared_);
        assert(isfinite(R2(k)) && R2(k)>=0 && R2(k)<=1+1e-8, ...
            'R2 de FOOOF fuera de rango.');
        pp=np.asarray(fm.peak_params_);
        if int64(py.len(pp)) == 0
            % Ajuste aperiódico valido SIN picos periodicos.
            pk=zeros(0,3);
        else
            vals=double(py.array.array('d',py.numpy.nditer(pp)));
            assert(mod(numel(vals),3)==0,'Formato inesperado de parametros de picos.');
            pk=reshape(vals,3,[])';
        end
        assert(size(pk,2)==3 && all(isfinite(pk(:))), ...
            'Parametros de picos mal formados o no finitos.');
        NumeroPicos(k)=size(pk,1);
        Picos{k}=pk;
        Estado(k)="AJUSTADO";
        if HACER_FIGURAS && numel(Ejemplos)<N_EJEMPLOS
            z.Sujeto=Sujeto(k); z.Electrodo=ElectrodoNombre(k);z.Momento=NombreMomento(k);
            z.Frec=FREQ_FIT; z.Original=yy(:);z.Full=yh(:);
            z.AP=double(py.array.array('d',py.getattr(fm,'_ap_fit'))); z.AP=z.AP(:);
            Ejemplos(end+1)=z; %#ok<SAGROW>
        end
    catch ME
        Offset(k)=nan;Exponent(k)=nan;R2(k)=nan;R2_Corr(k)=nan;
        ErrorFOOOF(k)=nan;NumeroPicos(k)=nan;Picos{k}=[];
        Estado(k)="ERROR_FOOOF";MensajeError(k)=string(ME.message);
    end
    if mod(k,500)==0, fprintf('FOOOF: %d/%d procesados\n',k,N);end
end

%% ===== VARIABLES POR BANDA: PICOS REALES DE FOOOF =====
% Bandas (limites continuos sin solapamiento, respetando la malla de 0.5 Hz):
% Delta parcial [2,4), Theta [4,8), Alfa [8,13), Beta [13,30), Gamma [30,40].
% El ajuste comienza en 2 Hz: Delta es PARCIAL (2-4 Hz nominales).
% FOOOF devuelve por pico [CF_Hz, PW_log10, BW_Hz].
% Cuando no se detecta pico: NPeaks=0, Peak=0, CF/BW=NaN.
% Cuando no se pudo ajustar el espectro: TODAS las variables son NaN.
bandNames = {'Delta','Theta','Alfa','Beta','Gamma'};
bandEdges = [2 4 8 13 30 40];
NPeaksBand = nan(N,5); PeakBand = nan(N,5);
CFBand = nan(N,5); BWBand = nan(N,5);
for k=1:N
    if Estado(k) ~= "AJUSTADO", continue; end
    pk=Picos{k};
    for b=1:5
        if b<5
            sel=pk(:,1)>=bandEdges(b) & pk(:,1)<bandEdges(b+1);
        else
            sel=pk(:,1)>=bandEdges(b) & pk(:,1)<=bandEdges(b+1);
        end
        subset=pk(sel,:);
        NPeaksBand(k,b)=size(subset,1);
        if isempty(subset)
            PeakBand(k,b)=0;
        else
            [~,im]=max(subset(:,2));
            PeakBand(k,b)=subset(im,2);
            CFBand(k,b)=subset(im,1);
            BWBand(k,b)=subset(im,3);
        end
    end
end
assert(all(sum(NPeaksBand(Estado=="AJUSTADO",:),2)==NumeroPicos(Estado=="AJUSTADO")), ...
    'Error de conteo: las bandas no reproducen todos los picos');

%% ================= TRES EXCEL INDEPENDIENTES =================
% Salida principal: SIEMPRE 16359 filas, sin eliminar ni reordenar momentos.
% Control final de forma ANTES de crear las tablas y los Excel.
assert(iscolumn(Sujeto) && iscolumn(Grupo) && iscolumn(Electrodo) && ...
    iscolumn(ElectrodoNombre) && iscolumn(Momento) && ...
    iscolumn(NombreMomento) && iscolumn(Fase) && ...
    iscolumn(Offset) && iscolumn(Exponent) && iscolumn(R2) && ...
    iscolumn(R2_Corr) && iscolumn(ErrorFOOOF) && ...
    iscolumn(NumeroPicos) && iscolumn(Estado) && iscolumn(NValidas), ...
    'Se han detectado vectores fila antes de construir RESULTADOS.');
Resultados=table(Sujeto,Grupo,Electrodo,ElectrodoNombre,Momento,NombreMomento,Fase, ...
    Offset,Exponent,R2,R2_Corr,ErrorFOOOF,NumeroPicos,Estado,NValidas, ...
    'VariableNames',{'Sujeto','Grupo','Electrodo_Num','Electrodo','Momento','EtiquetaMomento', ...
    'Fase','Offset','Exponent','R2','R2_Corr','Error','NumPicos','Estado','FrecuenciasValidas_2_40'});
% Mantener los nombres historicos Peak_Delta, Peak_Theta, etc.
% Cada fila corresponde a UNA combinacion sujeto-electrodo-momento.
for b=1:5
    name=bandNames{b};
    Resultados.(['Peak_' name])=PeakBand(:,b);
    Resultados.(['NPeaks_' name])=NPeaksBand(:,b);
    Resultados.(['CF_' name])=CFBand(:,b);
    Resultados.(['BW_' name])=BWBand(:,b);
end
assert(height(Resultados)==N,'Numero de filas de resultado incorrecto');
assert(all(cellfun(@(c) size(Resultados.(c),1)==N, Resultados.Properties.VariableNames)), ...
    'Dimensiones incorrectas en una columna de RESULTADOS.');
assert(numel(unique(strcat(Resultados.Sujeto, "_",string(Resultados.Electrodo_Num), "_",string(Resultados.Momento))))==N, ...
    'Claves repetidas en resultados');

% Salida QC: misma clave que resultados, causas y errores.
QC=table(Sujeto,Grupo,Electrodo,ElectrodoNombre,Momento,NombreMomento,Fase, ...
    NValidas,NNaN_Rango,NNaN_Global,EsCompleto,Estado,R2,R2_Corr,ErrorFOOOF,MensajeError, ...
    'VariableNames',{'Sujeto','Grupo','Electrodo_Num','Electrodo','Momento', ...
    'EtiquetaMomento','Fase','NValidas_2_40','NaN_2_40', ...
    'NaN_0_5_40','EspectroCompleto','Estado','R2','R2_Corr','Error','DetalleError'});
assert(sum(Estado=="AJUSTADO")+sum(Estado=="ERROR_FOOOF")+sum(~EsCompleto)==N);

% Desgloses por sujeto / electrodo / fase: recuentos, denominadores y porcentajes.
SU=unique(Sujeto,'sorted');EL=(1:19)';FA=["PRE_OA";"tSMS";"POST"];
QCSuj=tablaQC(Sujeto,SU,Estado);QCElec=tablaQC(Electrodo,EL,Estado);
QCFase=tablaQC(Fase,FA,Estado);
nombresElectrodosQC=string(canales(QCElec.Identificador));
QCElec.Electrodo=nombresElectrodosQC(:); % obligatorio: 19 x 1
assert(height(QCElec)==19 && numel(QCElec.Electrodo)==19, ...
    'Nombres de electrodos de auditoria mal dimensionados.');

% Tabla de TODOS los picos: archivo independiente, UNA hoja.
ppS=strings(0,1); ppG=zeros(0,1); ppE=zeros(0,1); ppEN=strings(0,1);
ppM=zeros(0,1);ppEtiqueta=strings(0,1);ppFase=strings(0,1);
ppId=zeros(0,1);ppBanda=strings(0,1);
ppCentro=zeros(0,1);ppAltura=zeros(0,1);ppAnchura=zeros(0,1);
for k=1:N
    pk=Picos{k}; if isempty(pk), continue; end
    n=size(pk,1);
    ppS=[ppS;repmat(Sujeto(k),n,1)]; %#ok<AGROW>
    ppG=[ppG;repmat(Grupo(k),n,1)]; %#ok<AGROW>
    ppE=[ppE;repmat(Electrodo(k),n,1)]; %#ok<AGROW>
    ppEN=[ppEN;repmat(ElectrodoNombre(k),n,1)]; %#ok<AGROW>
    ppM=[ppM;repmat(Momento(k),n,1)]; %#ok<AGROW>
    ppEtiqueta=[ppEtiqueta;repmat(NombreMomento(k),n,1)]; %#ok<AGROW>
    ppFase=[ppFase;repmat(Fase(k),n,1)]; %#ok<AGROW>
    ppId=[ppId;(1:n)']; %#ok<AGROW>
    ppCentro=[ppCentro;pk(:,1)]; %#ok<AGROW>
    ppAltura=[ppAltura;pk(:,2)]; %#ok<AGROW>
    ppAnchura=[ppAnchura;pk(:,3)]; %#ok<AGROW>
    for j=1:n
        cf=pk(j,1); band="Fuera_de_bandas";
        for b=1:5
            if (b<5 && cf>=bandEdges(b) && cf<bandEdges(b+1)) || ...
               (b==5 && cf>=bandEdges(b) && cf<=bandEdges(b+1))
                band=string(bandNames{b});break;
            end
        end
        ppBanda(end+1,1)=band; %#ok<AGROW>
    end
end
DetallePicos=table(ppS,ppG,ppE,ppEN,ppM,ppEtiqueta,ppFase,ppId,ppBanda,ppCentro,ppAltura,ppAnchura, ...
    'VariableNames',{'Sujeto','Grupo','Electrodo_Num','Electrodo','Momento','EtiquetaMomento', ...
    'Fase','IndicePico','Banda','CF_Hz','PW_log10','BW_Hz'});

nAjustados=nnz(Estado=="AJUSTADO");nIncompletos=nnz(Estado=="ESPECTRO_INCOMPLETO");
nSinPSD=nnz(Estado=="SIN_PSD_EN_RANGO");nErrores=nnz(Estado=="ERROR_FOOOF");
Resumen=table(N,nAjustados,nIncompletos,nSinPSD,nErrores,77,RANGO(1),RANGO(2), ...
    'VariableNames',{'EspectrosEsperados','Ajustados','Incompletos', ...
    'SinPSDEnRango','ErrorFOOOF','FrecuenciasRequeridas', ...
    'FrecuenciaMinHz','FrecuenciaMaxHz'});
Parametros=table(["FOOOF_python_version";"rango_Hz";"modo_aperiodico"; ...
    "peak_width_limits_Hz";"max_n_peaks";"min_peak_height";"peak_threshold"; ...
    "cobertura";"regla_R2";"potencia_entrada";"mapeo_temporal";"bandas_picos";"regla_pico_dominante";"cero_y_nan"], ...
    [string(versionFooof);string(sprintf("%.1f-%.1f",RANGO));string(MODO_AP); ...
     string(sprintf("%g-%g",ANCHO_PICOS));string(MAX_PICOS); ...
     string(MIN_ALTURA);string(UMBRAL_PICOS); ...
    "77/77 frecuencias";"sin exclusion por R2";"PSD lineal ORIGINAL sin IQR"; ...
    "1=PRE_OA,2:16=tSMS_1:15,17:21=POST_1:5"; ...
    "Delta_parcial=[2,4),Theta=[4,8),Alfa=[8,13),Beta=[13,30),Gamma=[30,40]"; ...
    "Maxima PW por banda; CF/BW del mismo pico. Peak_Delta solo 2-4 Hz"; ...
    "Sin pico: Peak=0,NPeaks=0,CF/BW=NaN; sin ajuste: todo NaN"], ...
    'VariableNames',{'Parametro','Valor'});

% Proteccion contra sobrescritura comprobada antes del bucle de ajustes.
assert(height(DetallePicos)==sum(NumeroPicos(Estado=="AJUSTADO")), ...
    'El detalle de picos no coincide con el recuento total.');
writetable(Resultados,OUT_RESULTADOS,'Sheet','RESULTADOS');
writetable(DetallePicos,OUT_PICOS,'Sheet','TODOS_LOS_PICOS');
writetable(Resumen,OUT_AUDITORIA,'Sheet','RESUMEN');
writetable(QC,OUT_AUDITORIA,'Sheet','ESTADO_POR_ESPECTRO');
writetable(QCSuj,OUT_AUDITORIA,'Sheet','POR_SUJETO');
writetable(QCElec,OUT_AUDITORIA,'Sheet','POR_ELECTRODO');
writetable(QCFase,OUT_AUDITORIA,'Sheet','POR_FASE');
QCGru=tablaQC(Grupo,[1;2],Estado);
writetable(QCGru,OUT_AUDITORIA,'Sheet','POR_GRUPO');
writetable(Parametros,OUT_AUDITORIA,'Sheet','PARAMETROS');


%% ================= FIGURAS QC =================
if HACER_FIGURAS
    r=R2(Estado=="AJUSTADO");
    if ~isempty(r)
        fig=figure('Visible','off','Color','w','Position',[100 100 1050 580]);
        histogram(r,35);xlabel('R^2');ylabel('Numero de espectros');
        title(sprintf('Calidad de ajuste FOOOF 2-40 Hz - n=%d',numel(r)));
        grid on;exportgraphics(fig,fullfile(CARPETA_FIG,'01_Histograma_R2.png'),'Resolution',300);close(fig);
    end
    fig=figure('Visible','off','Color','w','Position',[100 100 1050 580]);
    histogram(Exponent(Estado=="AJUSTADO"),40);xlabel('Exponente aperiodico');ylabel('Numero de espectros');
    title('Distribucion de exponentes');grid on;
    exportgraphics(fig,fullfile(CARPETA_FIG,'02_Distribucion_Exponentes.png'),'Resolution',300);close(fig);
    fig=figure('Visible','off','Color','w','Position',[100 100 1050 580]);
    histogram(ErrorFOOOF(Estado=="AJUSTADO"),40);xlabel('Error FOOOF');ylabel('Numero de espectros');
    title('Distribucion de errores');grid on;
    exportgraphics(fig,fullfile(CARPETA_FIG,'03_Distribucion_Errores.png'),'Resolution',300);close(fig);
    if height(DetallePicos)>0
        fig=figure('Visible','off','Color','w','Position',[100 100 1050 580]);
        histogram(DetallePicos.CF_Hz,38);xlabel('Frecuencia central real del pico (Hz)');ylabel('Numero de picos');
        title('Distribucion de centros de picos detectados por FOOOF');grid on;
        exportgraphics(fig,fullfile(CARPETA_FIG,'04_CentrosPicos.png'),'Resolution',300);close(fig);
    end
    for q=1:numel(Ejemplos)
        zz=Ejemplos(q);
        fig=figure('Visible','off','Color','w','Position',[100 100 1100 650]);
        plot(zz.Frec,zz.Original,'Color',[.4 .4 .4],'LineWidth',1.1);hold on;
        plot(zz.Frec,zz.Full,'Color',[.65 .15 .55],'LineWidth',2);
        plot(zz.Frec,zz.AP,'b--','LineWidth',1.5);
        xlabel('Frecuencia (Hz)');ylabel('log_{10}(PSD [\muV^2/Hz])');
        title(sprintf('%s - %s - %s',zz.Sujeto,zz.Electrodo,zz.Momento));
        legend({'Espectro original','Modelo completo','Componente aperiodico'},'Location','best');
        grid on;exportgraphics(fig,fullfile(CARPETA_FIG,sprintf('05_Ejemplo_Ajuste_%02d.png',q)),'Resolution',300);close(fig);
    end
end

%% ================= FIGURAS POR GRUPO Y TOTAL =================
% Estas son DESCRIPTIVAS, no contrastes estadisticos.
% Para evitar pseudorreplicacion en el grand average espectral:
% 1) media PSD lineal entre electrodos/momentos de cada sujeto (disponibles),
% 2) media de sujetos con IGUAL peso; se ajusta FOOOF de nuevo a cada media.
% No se utilizan picos a frecuencias prefijadas ni ruido aleatorio.
if HACER_FIGURAS
    etiquetas = {'Grupo 1','Grupo 2','Muestra total'};
    paleta = [1.00 .55 .41; .25 .82 .76; .57 .28 .68];
    fig=figure('Visible','off','Color','w','Position',[50 80 1600 540]);
    gruposSujeto=nan(numel(subjects),1);
    for is=1:numel(subjects)
        gruposSujeto(is)=Grupo(find(Sujeto==subjects(is),1));
    end
    for ig=1:3
        if ig<3
            sujetosGrupo=subjects(gruposSujeto==ig);
        else
            sujetosGrupo=subjects;
        end
        psdsSuj=nan(77,numel(sujetosGrupo));
        for is=1:numel(sujetosGrupo)
            selec=Sujeto==sujetosGrupo(is);
            psdsSuj(:,is)=mean(PSDfit(:,selec),2,'omitnan');
        end
        psdGrupo=mean(psdsSuj,2,'omitnan');
        assert(all(isfinite(psdGrupo) & psdGrupo>0), ...
            'PSD media incompleta en %s: comprobar cobertura por sujeto/frecuencia.',etiquetas{ig});
        fmMedia=fooofpy.FOOOF(py.list(ANCHO_PICOS),MAX_PICOS,MIN_ALTURA, ...
            UMBRAL_PICOS,MODO_AP,false);
        fmMedia.fit(np.array(FREQ_FIT'),np.array(psdGrupo'),py.list(RANGO));
        curvaOriginal=double(py.array.array('d',fmMedia.power_spectrum));
        curvaCompleta=double(py.array.array('d',fmMedia.fooofed_spectrum_));
        curvaAperiodica=double(py.array.array('d',py.getattr(fmMedia,'_ap_fit')));
        assert(numel(curvaOriginal)==77 && numel(curvaCompleta)==77 && ...
            numel(curvaAperiodica)==77 && all(isfinite(curvaCompleta)), ...
            'Ajuste del grand average incompleto en %s.',etiquetas{ig});
        subplot(1,3,ig);hold on;
        plot(FREQ_FIT,curvaOriginal(:),'Color',[.68 .68 .68],'LineWidth',1.2);
        plot(FREQ_FIT,curvaAperiodica(:),'b--','LineWidth',1.5);
        plot(FREQ_FIT,curvaCompleta(:),'Color',paleta(ig,:),'LineWidth',2.4);
        title(sprintf('%s (n=%d sujetos)',etiquetas{ig},numel(sujetosGrupo)));
        xlabel('Frecuencia (Hz)');ylabel('log_{10}(PSD [\muV^2/Hz])');
        xlim(RANGO);grid on;
        if ig==3
            legend({'PSD media observada','Componente aperiodica', ...
                'Modelo FOOOF de PSD media'},'Location','best');
        end
    end
    exportgraphics(fig,fullfile(CARPETA_FIG,'06_Grand_Average_FOOOF_Grupos_Total.png'), ...
        'Resolution',300);close(fig);

    % Tres informes descriptivos; n espectros != n sujetos.
    for ig=1:3
        if ig<3, sel=Estado=="AJUSTADO" & Grupo==ig;
        else, sel=Estado=="AJUSTADO"; end
        if ~any(sel),continue;end
        color=paleta(ig,:);
        if ig==1, suf='Grupo1';elseif ig==2,suf='Grupo2';else,suf='Total';end
        fig=figure('Visible','off','Color','w','Position',[80 70 1200 850]);
        subplot(2,2,1);
        histogram(Exponent(sel),40,'FaceColor',color);
        xlabel('Exponente aperiodico');ylabel('Numero de espectros');grid on;
        title('Distribucion de exponentes');
        subplot(2,2,2);
        histogram(R2(sel),30,'FaceColor',color);
        xlabel('R^2');ylabel('Numero de espectros');xlim([0 1]);grid on;
        title('Calidad del ajuste R^2');
        subplot(2,2,3);
        histogram(ErrorFOOOF(sel),35,'FaceColor',color);
        xlabel('Error FOOOF');ylabel('Numero de espectros');grid on;
        title('Error del modelo');
        subplot(2,2,4);
        if ig<3, ps=DetallePicos.Grupo==ig;
        else,ps=true(height(DetallePicos),1);end
        if any(ps)
            histogram(DetallePicos.CF_Hz(ps),2:1:40,'FaceColor',color);
        end
        xlim([2 40]);xlabel('Frecuencia central real (Hz)');
        ylabel('Numero de picos');grid on;
        title('Centros reales de picos detectados');
        sgtitle(sprintf('FOOOF - %s: %d espectros, %d sujetos', ...
            etiquetas{ig},nnz(sel),numel(unique(Sujeto(sel)))));
        exportgraphics(fig,fullfile(CARPETA_FIG, ...
            ['07_Reporte_Ajuste_' suf '.png']),'Resolution',300);close(fig);
    end
end

fprintf('\nFIN: Ajustados %d de %d. Incompletos %d, sin PSD %d, errores %d.\n', ...
    nAjustados,N,nIncompletos,nSinPSD,nErrores);
fprintf('Resultados: %s\nTodos los picos: %s\nAuditoria: %s\n',OUT_RESULTADOS,OUT_PICOS,OUT_AUDITORIA);

%% ================= FUNCIONES LOCALES =================
function Q = tablaQC(clave,ids,estado)
    ids=ids(:); n=numel(ids);total=zeros(n,1);ajustado=total;incompleto=total;sinPSD=total;errorFit=total;
    for j=1:n
        mask=(clave==ids(j));
        total(j)=nnz(mask);
        ajustado(j)=nnz(mask & estado=="AJUSTADO");
        incompleto(j)=nnz(mask & estado=="ESPECTRO_INCOMPLETO");
        sinPSD(j)=nnz(mask & estado=="SIN_PSD_EN_RANGO");
        errorFit(j)=nnz(mask & estado=="ERROR_FOOOF");
    end
    assert(all(total==ajustado+incompleto+sinPSD+errorFit));
    Q=table(ids,total,ajustado,incompleto,sinPSD,errorFit,100*(total-ajustado)./total, ...
        'VariableNames',{'Identificador','Total','Ajustados','Incompletos','SinPSD','ErrorFOOOF','PctNoAjustados'});
end