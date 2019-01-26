%% generate multiple datasets
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ % 
% specify starting number, ending number and file path
% assuming in path, you have four separate folders named: 
% saturation, pressure, production, rock

% *current configuration generates random layers only, well location is
% unaffected
PATH = '~/Desktop/caam-495-fluid-flow/New_Data/';
START = 1;
END = 2;
rng('shuffle') % randomize input based on clock

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ % 
% Solve layers of SPE10 with sequential and fully implicit solver
% This is a short example demonstrating how to solve any number of layers
% of the SPE10, model 2 problem using both a sequential and a fully
% implicit solver.
mrstModule add ad-core ad-blackoil spe10 blackoil-sequential mrst-gui

% Set up pressure and transport linear solvers
if ~isempty(mrstPath('agmg'))
    mrstModule add agmg
    psolver = AGMGSolverAD();
else
    psolver = BackslashSolverAD();
end
tsolver = GMRES_ILUSolverAD();
% Select layer 1
layers = 1; % change layer number
mrstModule add ad-core ad-blackoil blackoil-sequential spe10

% The base case for the model is 2000 days. This can be reduced to make the
% problem faster to run.
T = 2000*day;
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ % 


for NODE = START:END
    fprintf('Running %i/%i simulations...\n',NODE-START+1,END-START+1)
    fn_saturation = [PATH,'saturation/saturaton',num2str(NODE),'.txt'];
    fn_pressure = [PATH,'pressure/pressure',num2str(NODE),'.txt'];
    fn_production = [PATH,'production/production',num2str(NODE),'.txt'];
    fn_rock = [PATH,'rock/rock',num2str(NODE),'.txt'];
    
% hyperparameter inputs %
% multipliers: 
%[ layer porosity, layer permeability, well rate, producer rate]
multiplier = [1.0, 1.0, 1.0, 1.0]; % randomize multiplier
% well locations (fix for now)
well_loc.w1 = [1,1];    
well_loc.w2 = [60,1];   
well_loc.w3 = [60,220];
well_loc.w4 = [1,220];  
% producer locations (fix for now)
prod_loc    = [30,110];
% mixing layers [layer1, layer2, layer1 percent, layer2 percent]
a = rand;
mixed_layers = [randi(50)+35,randi(50)+35,a,1-a];
fprintf('Selecting layers %i and %i...\n', mixed_layers(1),mixed_layers(2))
[state, model, schedule] = setupSPE10_AD('layers', layers, 'dt', 30*day, ...
                'T',  T, 'well_loc',well_loc, 'prod_loc', prod_loc,...
                'multiplier',multiplier,'mixed_layers',mixed_layers);

% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ % 


%% Sequential
% Set up the sequential model
seqModel = getSequentialModelFromFI(model, 'pressureLinearSolver', psolver,....
                                           'transportLinearSolver', tsolver);
% We set up a timestep selector that aims for timesteps where the
% maximum saturation change is equal to a fixed value.
stepSel = StateChangeTimeStepSelector('targetProps', {'s'},...
                                      'targetChangeAbs', 0.25);
% Run problem
solver = NonLinearSolver('timeStepSelector', stepSel);
[wsSeq, statesSeq, repSeq] = simulateScheduleAD(state, seqModel, schedule, 'NonLinearSolver', solver);

%% Simulate fully implicit 
% Solve the fully implicit version of the problem, with a CPR
% preconditioner that uses the same linear solver for the pressure as the
% sequential solver.
% solver.timeStepSelector.reset();
% solver.LinearSolver = CPRSolverAD('ellipticSolver', psolver);

%[wsFIMP, statesFIMP, repFIMP] = simulateScheduleAD(state, model, schedule, 'NonLinearSolver', solver);

%% Plot simulation time taken
% The sequential solver can be much faster than the fully implicit solver
% for certain problems.
% figure(1); clf; hold on
% plot(cumsum(repSeq.SimulationTime)/60, 'b-*')
% % plot(cumsum(repFIMP.SimulationTime)/60, 'r-o')
% ylabel('Simulation time [minutes]')
% xlabel('Control step #')
% % legend('Sequential implicit', 'Fully implicit')
% legend('Sequential implicit')

%% Plot the results in interactive viewers Sequential
% G = model.G;
% W = schedule.control(1).W;
% 
% % Plot the well curves
% plotWellSols({wsSeq}, cumsum(schedule.step.val), ...
%     'datasetnames', {'Sequential'}, 'field', 'qOs')
% 
% % Plot reservoir quantities
% figure;
% plotToolbar(G, statesSeq);
% axis equal tight
% view(90, 90);
% plotWell(G, W);
% title('Sequential implicit')


%% Fully-implicit Plots
% figure;
% plotToolbar(G, statesFIMP);
% axis equal tight
% view(90, 90);
% plotWell(G, W);
% title('Fully-implicit')

% plotWellSols({wsSeq, wsFIMP}, cumsum(schedule.step.val), ...
%             'datasetnames', {'Sequential', 'FIMP'}, 'field', 'qOs')
%% 
% figure;
% subplot(1, 2, 1)
% bar([sum(repFIMP.SimulationTime), sum(repSeq.SimulationTime)]/60)
% set(gca, 'XTickLabel', {'Fully-implicit', 'Sequential'})
% ylabel('Simulation time [minutes]')
% set(gca, 'FontSize', 18)
% axis tight
% subplot(1, 2, 2)
% plotCellData(model.G, statesSeq{60}.s(:, 1), 'edgecolor', 'none')
% title('S_w')
% set(gca, 'FontSize', 18)
% 
% axis equal tight off


%% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ % 
% Data pulling
% each sequence of data is stored in row format
% number of rows = number of sequences
fprintf('Saving data...')

tic
% store pressures
dlmwrite(fn_pressure,statesSeq{1}.pressure')
for i = 2:75
    dlmwrite(fn_pressure,statesSeq{i}.pressure','-append')
end
% store oil pressures
dlmwrite(fn_saturation,statesSeq{1}.s(:,1)')
for i = 2:75
    dlmwrite(fn_saturation,statesSeq{i}.s(:,1)','-append')
end

% rock contains first permeability and then porosity
dlmwrite(fn_rock,model.rock.poro')
for i = 1:3
    dlmwrite(fn_rock,model.rock.perm','-append')
end

% production(m^3/s) together with time stamps(s)
file = fopen(fn_production,'w');
schedules = cumsum(schedule.step.val);
for i = 1:75
    fprintf(file,'%e\t',schedules(i));
    for j = 1:4
        fprintf(file,'%e\t',-wsSeq{i}(j).qOr);
    end
    fprintf(file,'\n');
end
fclose(file);

toc

end
%% Copyright notice

% <html>
% <p><font size="-1">
% Copyright 2009-2018 SINTEF ICT, Applied Mathematics.
% </font></p>
% <p><font size="-1">
% This file is part of The MATLAB Reservoir Simulation Toolbox (MRST).
% </font></p>
% <p><font size="-1">
% MRST is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% </font></p>
% <p><font size="-1">
% MRST is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% </font></p>
% <p><font size="-1">
% You should have received a copy of the GNU General Public License
% along with MRST.  If not, see
% <a href="http://www.gnu.org/licenses/">http://www.gnu.org/licenses</a>.
% </font></p>
% </html>
