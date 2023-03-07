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
    beta1Min = 0.0001;
    beta1Max = 0.04;
    beta2Min = 0.000000001;
    beta2Max = 0.0001;
    simuls = 10000;
    varsNumber = 3;
    savePath = "/Volumes/Aortas/Aorta/results/mathmodel/model3/1CellCycle/SimulationP3_1WOLoss/";

    X = lhsdesign(simuls,varsNumber);

    for i = 1:simuls
        cycle1 = X(i,1)*(cycleMax - cycleMin) + cycleMin;
        beta1 = X(i,2)*(beta1Max - beta1Min) + beta1Min;
        %beta1 = 0;
        %beta2 = 0;
        beta2 = X(i,3)*(beta2Max - beta2Min) + beta2Min;
        [T,Y] = ode15s(@(T,Y)msc(T,Y,fitresultNorm, cycle1,beta1, beta2), ...
                [timeStart timeEnd], [cellsStart]); 
        filename =  savePath + "/result" + i + ".csv";
        betas1 = beta1*ones(length(T),1);
        betas2 = beta2*ones(length(T),1);
        cycle1s = cycle1*ones(length(T),1);
        writematrix([T,Y, cycle1s, betas1, betas2],filename);
    end
end


function dy=msc(t,y,fitresultNorm, cycle1, beta1, beta2)
    a = fitresultNorm(t);
    growth_rate = 24/cycle1;  
    dy=zeros(2,1); % initiate array
    
    dy = (((a* growth_rate)/(1+(beta2*y))) - beta1)* y;  % Cells
end