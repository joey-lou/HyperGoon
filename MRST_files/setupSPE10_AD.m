function [state, model, schedule]  = setupSPE10_AD(varargin)
    well_loc.w1 = [1,1];
    well_loc.w2 = [60,1];
    well_loc.w3 = [60,220];
    well_loc.w4 = [1,220];
    opt = struct('layers', 1:85, ...
                 'dt',      30*day, ...
                 'T',       2000*day, ...
                 'minporo', 0.01, ...
                 'multiplier', [1.0, 1.0, 1.0, 1.0],...
                 'prod_loc', [20,100],...
                 'well_loc', well_loc,...
                 'mixed_layers', [36,36,0.5,0.5]  );
    opt = merge_options(opt, varargin{:});
    
    mrstModule add spe10 ad-props ad-blackoil ad-core
    


    srw = 0.2;
    sro = 0.2;
    pRef = 6000*psia;
    % impose multiplier for all permeabilities
    % Fluid relative permeabilities (vary permeabilities)
    fluid.krW = coreyPhaseRelpermAD(2, srw, 1, srw + sro);
    fluid.krO = coreyPhaseRelpermAD(2, sro, 1, srw + sro);
    
    % Water props
    bW0 = 1./(1.01);
    fluid.bW = @(p) bW0*exp((p - pRef)*3e-6/psia);
    fluid.muW = @(p) 0*p + 0.3*centi*poise;
    
    % Oil props
    p = [300; 800; 8000; 8001]*psia;
    b = 1./[1.05; 1.02; 1.01; 1.01];
    mu = [2.85; 2.99; 3; 3]*centi*poise;
    [fluid.bO, fluid.muO] = tableByPressureLinearAD(p, b, mu);
    
    
    fluid.rhoWS = 64*pound/(ft^3);
    fluid.rhoOS = 53*pound/(ft^3);

    % Rock compressibility
    cR = 1e-6/psia;
    fluid.cR = cR;
    fluid.pvMultR = @(p)(1 + cR.*(p-pRef));
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ % 
    % changing permeability and porosity
    rock_1 = getSPE10rock(opt.mixed_layers(1));
    rock_2 = getSPE10rock(opt.mixed_layers(2));
    rock.perm = opt.mixed_layers(3)*rock_1.perm+...
        opt.mixed_layers(4)*rock_2.perm;
    rock.poro = opt.mixed_layers(3)*rock_1.poro+...
        opt.mixed_layers(4)*rock_2.poro;
    % (1) scale rock porosity by a multiplier [+-10%]
    rock.poro = rock.poro*opt.multiplier(1);
    % (2) scale rock permeability by a multiplier [+-10%]
    rock.perm = rock.perm*opt.multiplier(2);
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ %     
    
    % Compute pore volume fraction of the full model
    volFrac = sum(rock.poro)/1.9141e+05;
    rock.poro(rock.poro < opt.minporo) = opt.minporo;
    
    % Grid
    cartDims = [  60,  220,   numel(opt.layers)];
    physDims = [1200, 2200, 2*cartDims(end)] .* ft();

    G = cartGrid(cartDims, physDims);
    try
        mrstModule add libgeometry
        G = mcomputeGeometry(G);
    catch
        G = computeGeometry(G);
    end
    model = TwoPhaseOilWaterModel(G, rock, fluid, 'gravity', [0, 0, 9.80665]);
    model.minimumPressure = 1000*psia;
    
    state = initResSol(G, pRef, [srw, 1-srw]);
    
    % Wells
    % (3) impose [+-10%] 
    makeProd = @(W, name, I, J) verticalWell(W, G, rock, I, J, [],...
        'Name', name, 'radius', 5*inch, 'sign', -1, 'Type', 'bhp',...
        'Val', 4000*psia*opt.multiplier(3), 'comp_i', [.5, .5]);
    W = [];
    % change well location
    W = makeProd(W, 'P1', opt.well_loc.w1(1), opt.well_loc.w1(2));
    W = makeProd(W, 'P2', opt.well_loc.w2(1), opt.well_loc.w2(2));
    W = makeProd(W, 'P3', opt.well_loc.w3(1), opt.well_loc.w3(2));
    W = makeProd(W, 'P4', opt.well_loc.w4(1), opt.well_loc.w4(2));
    % (4) impose [+-10%] 
    W = verticalWell(W, G, rock, opt.prod_loc(1), opt.prod_loc(2), [], 'Name', 'I1', 'radius', 5*inch, ...
        'Type', 'rate', 'Val', volFrac*5000*stb/day*opt.multiplier(4), 'comp_i', [1, 0], 'Sign', 1);

    dt = rampupTimesteps(opt.T, opt.dt);
    
    schedule = simpleSchedule(dt, 'W', W);
end
