% ---------------------------------------------------------
% ------------------ add pulseq objects -------------------
% ----------------- readout: spiral (SPI) -----------------
% ---------------------------------------------------------

% ----- sequence loop counter: -----
% loop_rf_inc -> rf spoiling
% loop_kz -> z partitions
% loop_NR -> repetitions

%% init loop counters
if ~exist('loop_rf_inc', 'var')
    loop_rf_inc = 0;
end
if ~exist('loop_kz', 'var')
    loop_kz = 1;
end
if loop_NR<1
    temp_loop = 1;
else
    temp_loop = loop_NR;
end

%% rf spoilng
if strcmp(SPI.spoil_rf_mode, 'import')
    temp_phase  = SPI.mrf_import.rf_phases(temp_loop);
else
    temp_phase  = SPI.spoil_rf_inc * loop_rf_inc ^ SPI.spoil_rf_pow;
    loop_rf_inc = loop_rf_inc + 1;
end
SPI.rf(temp_loop).phaseOffset = mod(temp_phase,        2*pi);
SPI.adc.phaseOffset           = mod(temp_phase + pi/2, 2*pi); % shift with pi/2: Mx -> real axis
clear temp_phase;

%% segment label extension for GE scanners
if flag_GE==1
    if loop_NR<1
        seq.addBlock(mr.makeLabel('SET', 'TRID', 0));
    else
        seq.addBlock(mr.makeLabel('SET', 'TRID', 1));
    end
end

%% LIN label extension for United Imaging scanners
if flag_UI==1
    if loop_NR>0
        if ~exist('loop_lin', 'var')
            loop_lin = 1;
        end
        seq.addBlock(mr.makeLabel('SET','LIN', loop_lin));
        loop_lin = loop_lin + 1;
    end
end

%% add sequence blocks
seq.addBlock(SPI.rf(temp_loop), SPI.gz);
seq.addBlock(SPI.gz_reph(loop_kz));
seq.addBlock(SPI.TE_delay(temp_loop));
if loop_NR<1
    seq.addBlock(SPI.gx(SPI.phi_id(temp_loop)), SPI.gy(SPI.phi_id(temp_loop)), SPI.gz_spoil(loop_kz));           % dummy loop 
else
    seq.addBlock(SPI.gx(SPI.phi_id(temp_loop)), SPI.gy(SPI.phi_id(temp_loop)), SPI.gz_spoil(loop_kz), SPI.adc);  % adc loop
end
seq.addBlock(SPI.TR_delay(temp_loop));
clear temp_loop