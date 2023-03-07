%number of total cells

function aortic_TotalCells
    rng('default') % for reproducibility of normrnd

    close all;
    % pull angle distributions, mitotic normalization, and cell length from
    % observed data:

    %fit func growth
    positive = readtable("/Volumes/Aortas/Aorta/results/cellSize/KI67ActiveP5_60.csv");
    [xData, yData] = prepareCurveData(positive.ages, positive.active);
    ft = 'pchipinterp';
    [fitresultNorm, gof] = fit( xData, yData, ft, 'Normalize', 'on' );

    cellsStart = 20608; %number of cells
    timeStart = 3;
    timeEnd = 60;
    cycleMin = 5;
    cycleMax = 80;
    betaMin = 0.0001;
    betaMax = 0.04;
    simuls = 1000;
    varsNumber = 1;
    savePath = "/Volumes/Aortas/Aorta/results/mathmodel/model1/1CellCycle/SimulationP3_1WOLoss/";


    X = lhsdesign(simuls,varsNumber);

    for i = 1:simuls
        cycle1 = X(i,1)*(cycleMax - cycleMin) + cycleMin;
        %beta = X(i,2)*(betaMax - betaMin) + betaMin;
        beta = 0
        [T,Y] = ode15s(@(T,Y)msc(T,Y,fitresultNorm, cycle1,beta), ...
                [timeStart timeEnd], [cellsStart]); 
        filename =  savePath + "/result" + i + ".csv";
        betas = beta*ones(length(T),1);
        cycle1s = cycle1*ones(length(T),1);
        writematrix([T,Y, cycle1s, betas],filename);
    end
end


function dy=msc(t,y,fitresultNorm, cycle1,beta)
    a = fitresultNorm(t);
    growth_rate = 24/cycle1;
    dy=zeros(2,1); % initiate array

    dy = ((a* growth_rate) - beta) * y;  % Cells
end