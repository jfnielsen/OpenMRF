% ---------------------------------------------------------
% ------------------ add pulseq objects -------------------
% ----------------- readout: spiral (SPI) -----------------
% -------------------- Noise pre-scans --------------------
% ---------------------------------------------------------

if isfield(SPI, 'Nnoise')
    for loop_noise = 1:SPI.Nnoise
        if flag_GE==1
            seq.addBlock(mr.makeLabel('SET', 'TRID', 1));
        end
        seq.addBlock(SPI.adc, ceil((mr.calcDuration(SPI.adc) + 1e-3) / system.blockDurationRaster) * system.blockDurationRaster);
    end
end
