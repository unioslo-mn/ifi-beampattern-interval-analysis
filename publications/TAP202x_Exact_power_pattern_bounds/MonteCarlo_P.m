function value = MonteCarlo_P(bp,options)

    arguments
        bp
        options.theta   (1,1)   {mustBeNumeric} = nan()
    end
    
    if isnan(options.theta)
        thetas = bp.BeamAngles;
    else
        thetas = options.theta;
    end
    
    % Calculate some needed quantities
    M = bp.Arrays.ElCount;
    w_apod = bp.Arrays.TaperWeights;
    k_abs = bp.Arrays.WaveNumber;
    lambda = bp.Arrays.WaveLength;
    elRad = bp.Arrays.ElDiameter/2;
    taperAngs = bp.Arrays.TaperAngles;
    k_s = bp.Arrays.SteeringVector;
    rx_0 = bp.Arrays.ElPosX;
    ry_0 = bp.Arrays.ElPosY;
    
    % Get Monte-Carlo structure
    g_I = bp.Arrays.GainInterval;
    ph_I = bp.Arrays.PhaseInterval;
    rx_I = bp.Arrays.PosXInterval;
    ry_I = bp.Arrays.PosYInterval;
    psi_I = bp.Arrays.OrientInterval;
    g_MC = g_I.Infimum + g_I.Width .* rand(M, 1);
    phi_MC = ph_I.Infimum + ph_I.Width .* rand(M, 1);
    rx_MC = rx_I.Infimum + rx_I.Width .* rand(M, 1);
    ry_MC = ry_I.Infimum + ry_I.Width .* rand(M, 1);
    psi_MC = psi_I.Infimum + psi_I.Width .* rand(M, 1);
    
    rand_C_prime = rand(M, M) .* exp( 1j* 2*pi * rand(M, M) );
    rand_C_prime( 1 : M + 1 : M^2) = 1;
    C_MC = bp.Arrays.CouplingAmp .* rand_C_prime;

   
    value = zeros(length(thetas), 1);
    for t = 1 : length(thetas) % Loop over angles
        k = k_abs * ([sin(thetas(t)), cos(thetas(t))] );
            
        sum_m = 0;
        for m = 1:M
            sum_m = sum_m + w_apod(m) * C_MC(m,:).' * exp( -1j*( k_s(1)*rx_0(m) + k_s(2)*ry_0(m) ));
        end

        delta_Psi = thetas(t) + (-psi_MC);
        Ec_ampl = g_MC .* bp.Arrays.getElDirectivity(delta_Psi);

        phase = sum(k .* [rx_MC, ry_MC],2) + phi_MC;
        Ec_ph = exp(1j.*phase);

        BP = sum( Ec_ampl .* Ec_ph .* sum_m ) ; 

        
        value(t) = abs(BP).^2;
    end    

end

