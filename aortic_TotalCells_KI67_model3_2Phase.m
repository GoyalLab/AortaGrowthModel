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
    cycle1Min = 5;
    cycle1Max = 30;
    cycle2Min = 20;
    cycle2Max = 80;
    beta1Min = 0.0001;
    beta1Max = 0.04;
    beta2Min = 0.000000001;
    beta2Max = 0.0001;
    simuls = 10000;
    varsNumber = 4;
    savePath = "/Volumes/Aortas/Aorta/results/mathmodel/model3/2CellCycle/SimulationP3_2/";

    X = lhsdesign(simuls,varsNumber);

    for i = 1:simuls
        cycle1 = X(i,1)*(cycle1Max - cycle1Min) + cycle1Min;
        cycle2 = X(i,2)*(cycle2Max - cycle2Min) + cycle2Min;
        beta1 = X(i,3)*(beta1Max - beta1Min) + beta1Min;
        beta2 = X(i,4)*(beta2Max - beta2Min) + beta2Min;
        [T,Y] = ode15s(@(T,Y)msc(T,Y,fitresultNorm, cycle1, cycle2, beta1, beta2), ...
                [timeStart timeEnd], [cellsStart]); 
        filename =  savePath + "/result" + i + ".csv";
        betas1 = beta1*ones(length(T),1);
        betas2 = beta2*ones(length(T),1);
        cycle1s = cycle1*ones(length(T),1);
        cycle2s = cycle2*ones(length(T),1);
        writematrix([T,Y, cycle1s, cycle2s, betas1, betas2],filename);
    end  
end


function dy=msc(t,y,fitresultNorm, cycle1, cycle2, beta1, beta2)  
    a = fitresultNorm(t);
    if t <= 14
        growth_rate = 24/cycle1;
    else
        growth_rate = 24/cycle2;
    end 
    dy=zeros(2,1); % initiate array
    
    dy = (((a* growth_rate)/(1+(beta2*y))) - beta1)* y;
end