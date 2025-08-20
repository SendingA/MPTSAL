function aoa_out = map_pdoa_to_aoa(pdoa_in)
    PDOA_LUT_1_CH9 = [...
        -3.3905, -1.5708;
        -3.2998, -1.4661;
        -3.1978, -1.3614;
        -3.0790, -1.2566;
        -2.9326, -1.1519;
        -2.7431, -1.0472;
        -2.5136, -0.9425;
        -2.2552, -0.8378;
        -1.9790, -0.7330;
        -1.6963, -0.6283;
        -1.4002, -0.5236;
        -1.1205, -0.4189;
        -0.8410, -0.3142;
        -0.5628, -0.2094;
        -0.2742, -0.1047;
         0.0000,  0.0000;
         0.2770,  0.1047;
         0.5553,  0.2094;
         0.8339,  0.3142;
         1.1050,  0.4189;
         1.3498,  0.5236;
         1.5912,  0.6283;
         1.8036,  0.7330;
         2.0040,  0.8378;
         2.1975,  0.9425;
         2.3780,  1.0472;
         2.5190,  1.1519;
         2.6505,  1.2566;
         2.7652,  1.3614;
         2.8451,  1.4661;
         2.8932,  1.5708];
    
    pdoa_list = PDOA_LUT_1_CH9 (:,1);
    aoa_list = PDOA_LUT_1_CH9 (:,2);
    if pdoa_in <= pdoa_list(1)
        aoa_out = aoa_list(1);
        return;
    elseif pdoa_in >= pdoa_list(end)
        aoa_out = aoa_list(end);
        return;
    end

    for i = 1:length(pdoa_list)-1
        if pdoa_in >= pdoa_list(i) && pdoa_in <= pdoa_list(i+1)
            pdoa0 = pdoa_list(i);
            pdoa1 = pdoa_list(i+1);
            aoa0 = aoa_list(i);
            aoa1 = aoa_list(i+1);
            aoa_out = aoa0 + (pdoa_in - pdoa0) * (aoa1 - aoa0) / (pdoa1 - pdoa0);
            return;
        end
    end

    error('PDOA value out of LUT range.');
end
