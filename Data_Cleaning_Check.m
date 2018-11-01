%%
clear all;close all;clc

speaker = 28;
for ar= 1:19;
    file = ['SCI' num2str(speaker) '_ar' num2str(ar)];
    DisplacementFigure = (['SCI' num2str(speaker) '/Displacement Figures/' file]);
    openfig(DisplacementFigure);
    
    MinCheckDorsum = .22;
    MaxCheckDorsum = .569;
    
    MinCheckRoot = .1;
    MaxCheckRoot = .574;
    
    MinCheckBlade = .1;
    MaxCheckBlade = .5936;
    
    
    DataPointsRef = findobj(subplot(2,2,1),'Type','line');
    xRef=get(DataPointsRef,'Xdata');
    yRef=get(DataPointsRef,'Ydata');
    
    numericCellsxRoot = xRef(1);
    numericVectorxRootRef = cell2mat(numericCellsxRoot);
    rootxRef = double(numericVectorxRootRef);
    
    numericCellsxDorsum = xRef(2);
    numericVectorxDorsumRef = cell2mat(numericCellsxDorsum);
    dorsumxRef = double(numericVectorxDorsumRef);
    
    numericCellsxBlade = xRef(3);
    numericVectorxBladeRef = cell2mat(numericCellsxBlade);
    bladexRef = double(numericVectorxBladeRef);
    
    numericCellsyRoot = yRef(1);
    numericVectoryRootRef = cell2mat(numericCellsyRoot);
    rootyRef = double(numericVectoryRootRef);
    
    numericCellsyDorsum = yRef(2);
    numericVectoryDorsumRef = cell2mat(numericCellsyDorsum);
    dorsumyRef = double(numericVectoryDorsumRef);
    
    numericCellsyBlade = yRef(3);
    numericVectoryBladeRef = cell2mat(numericCellsyBlade);
    bladeyRef = double(numericVectoryBladeRef);
    
    RefMaxRootx = max(rootxRef);
    RefMaxDorsumx = max(dorsumxRef);
    RefMaxBladex = max(bladexRef);
    
    RefMinRootx = min(rootxRef);
    RefMinDorsumx = min(dorsumxRef);
    RefMinBladex = min(bladexRef);
    
    RefMaxRooty = max(rootyRef);
    RefMaxDorsumy = max(dorsumyRef);
    RefMaxBladey = max(bladeyRef);
    
    RefMinRooty = min(rootyRef);
    RefMinDorsumy = min(dorsumyRef);
    RefMinBladey = min(bladeyRef);
    
    RootLength = RefMaxRootx - RefMinRootx;
    DorsumLength = RefMaxDorsumx -RefMinDorsumx;
    BladeLength = RefMaxBladex - RefMinBladex;
    
    RootHeight = RefMaxRooty - RefMinRooty;
    DorsumHeight = RefMaxDorsumy -RefMinDorsumy;
    BladeHeight = RefMaxBladey - RefMinBladey;
    %%%%%%%%%%%%%%%%%%%%%%%%
    DataPointsA = findobj(subplot(2,2,2),'Type','line');
    DataPointsR = findobj(subplot(2,2,3),'Type','line');
    
    xa=get(DataPointsA,'Xdata');
    ya=get(DataPointsA,'Ydata');
    
    xr=get(DataPointsR,'Xdata');
    yr=get(DataPointsR,'Ydata');
    
    numericCells = xa;
    numericVectorxa = cell2mat(numericCells);
    
    numericCells = ya;
    numericVectorya = cell2mat(numericCells);
    
    numericCells = xr;
    numericVectorxr = cell2mat(numericCells);
    
    numericCells = yr;
    numericVectoryr = cell2mat(numericCells);
    
    xatrans = transpose(numericVectorxa);
    rootxa = (xatrans(1:3,1));
    dorsumxa = (xatrans(1:3,2));
    bladexa = (xatrans(1:3,3));
    
    yatrans = transpose(numericVectorya);
    rootya = (yatrans(1:3,1));
    dorsumya = (yatrans(1:3,2));
    bladeya = (yatrans(1:3,3));
    
    xrtrans = transpose(numericVectorxr);
    rootxr = (xrtrans(1:3,1));
    dorsumxr = (xrtrans(1:3,2));
    bladexr = (xrtrans(1:3,3));
    
    yrtrans = transpose(numericVectoryr);
    rootyr = (yrtrans(1:3,1));
    dorsumyr = (yrtrans(1:3,2));
    bladeyr = (yrtrans(1:3,3));
    
    %%%%%%%%%%%%%%
    
    r1xaroot = rootxa(1);
    r1yaroot = rootya(1);
    r2xaroot = rootxa(2);
    r2yaroot = rootya(2);
    
    r1r2roota = [[r1xaroot, r1yaroot];[r2xaroot, r2yaroot]];
    
    r1r2DistRoota = pdist(r1r2roota, 'euclidean');
    
    r3xaroot = rootxa(3);
    r3yaroot = rootya(3);
    r2xaroot = rootxa(2);
    r2yaroot = rootya(2);
    
    r3r2roota = [[r3xaroot, r3yaroot];[r2xaroot, r2yaroot]];
    
    r3r2DistRoota = pdist(r3r2roota, 'euclidean');
    
    SegmentLengthRoota = abs (r1r2DistRoota + r3r2DistRoota);
    
    RootChecka = SegmentLengthRoota/RootLength;
    
    r1xadorsum = dorsumxa(1);
    r1yadorsum = dorsumya(1);
    r2xadorsum = dorsumxa(2);
    r2yadorsum = dorsumya(2);
    
    r1r2dorsuma = [[r1xadorsum, r1yadorsum];[r2xadorsum, r2yadorsum]];
    
    r1r2DistDorsuma = pdist(r1r2dorsuma, 'euclidean');
    
    r3xadorsum = dorsumxa(3);
    r3yadorsum = dorsumya(3);
    r2xadorsum = dorsumxa(2);
    r2yadorsum = dorsumya(2);
    
    r3r2dorsuma = [[r3xadorsum, r3yadorsum];[r2xadorsum, r2yadorsum]];
    
    r3r2DistDorsuma = pdist(r3r2dorsuma, 'euclidean');
    
    SegmentLengthDorsuma = abs (r1r2DistDorsuma + r3r2DistDorsuma);
    
    DorsumChecka = SegmentLengthDorsuma/DorsumLength;
    
    r1xablade = bladexa(1);
    r1yablade = bladeya(1);
    r2xablade = bladexa(2);
    r2yablade = bladeya(2);
    
    r1r2bladea = [[r1xablade, r1yablade];[r2xablade, r2yablade]];
    
    r1r2DistBladea = pdist(r1r2bladea, 'euclidean');
    
    r3xablade = bladexa(3);
    r3yablade = bladeya(3);
    r2xablade = bladexa(2);
    r2yablade = bladeya(2);
    
    r3r2Bladea = [[r3xablade, r3yablade];[r2xablade, r2yablade]];
    
    r3r2DistBladea = pdist(r3r2Bladea, 'euclidean');
    
    SegmentLengthBladea = abs (r1r2DistBladea + r3r2DistBladea);
    
    BladeChecka = SegmentLengthBladea/BladeLength;
    
    count = 0;
    
    if (MinCheckRoot < RootChecka) && (RootChecka < MaxCheckRoot)
        count = count + 1;
    else
        fprintf('BAD /a/ Frame Root Displacement | Production: %i\n',ar)
    end
    
    if (MinCheckDorsum < DorsumChecka) && ( DorsumChecka < MaxCheckDorsum)
        count = count + 1;
    else
        fprintf('BAD /a/ Frame Dorsum Displacements | Production: %i\n',ar)
    end
    
    if (MinCheckBlade < BladeChecka) && (BladeChecka < MaxCheckBlade)
        count = count + 1;
    else
        fprintf('BAD /a/ Frame Blade Displacements | Production: %i\n',ar)
    end
    
    if count == 3;
        fprintf('All /a/ Frame Displacements are GOOD\n')
    end
    %%%%%%%%%%%%%
    
    r1xrroot = rootxr(1);
    r1yrroot = rootyr(1);
    r2xrroot = rootxr(2);
    r2yrroot = rootyr(2);
    
    r1r2rootr = [[r1xrroot, r1yrroot];[r2xrroot, r2yrroot]];
    
    r1r2DistRootr = pdist(r1r2rootr, 'euclidean');
    
    r3xrroot = rootxr(3);
    r3yrroot = rootyr(3);
    r2xrroot = rootxr(2);
    r2yrroot = rootyr(2);
    
    r3r2rootr = [[r3xrroot, r3yrroot];[r2xrroot, r2yrroot]];
    
    r3r2DistRootr = pdist(r3r2rootr, 'euclidean');
    
    SegmentLengthRootr = abs (r1r2DistRootr + r3r2DistRootr);
    
    RootCheckr = SegmentLengthRootr/RootLength;
    
    r1xrdorsum = dorsumxr(1);
    r1yrdorsum = dorsumyr(1);
    r2xrdorsum = dorsumxr(2);
    r2yrdorsum = dorsumyr(2);
    
    r1r2dorsumr = [[r1xrdorsum, r1yrdorsum];[r2xrdorsum, r2yrdorsum]];
    
    r1r2DistDorsumr = pdist(r1r2dorsumr, 'euclidean');
    
    r3xrdorsum = dorsumxr(3);
    r3yrdorsum = dorsumyr(3);
    r2xrdorsum = dorsumxr(2);
    r2yrdorsum = dorsumyr(2);
    
    r3r2dorsumr = [[r3xrdorsum, r3yrdorsum];[r2xrdorsum, r2yrdorsum]];
    
    r3r2DistDorsumr = pdist(r3r2dorsumr, 'euclidean');
    
    SegmentLengthDorsumr = abs (r1r2DistDorsumr + r3r2DistDorsumr);
    
    DorsumCheckr = SegmentLengthDorsumr/DorsumLength;
    
    r1xrblade = bladexr(1);
    r1yrblade = bladeyr(1);
    r2xrblade = bladexr(2);
    r2yrblade = bladeyr(2);
    
    r1r2blader = [[r1xrblade, r1yrblade];[r2xrblade, r2yrblade]];
    
    r1r2DistBlader = pdist(r1r2blader, 'euclidean');
    
    r3xrblade = bladexr(3);
    r3yrblade = bladeyr(3);
    r2xrblade = bladexr(2);
    r2yrblade = bladeyr(2);
    
    r3r2Blader = [[r3xrblade, r3yrblade];[r2xrblade, r2yrblade]];
    
    r3r2DistBlader = pdist(r3r2Blader, 'euclidean');
    
    SegmentLengthBlader = abs (r1r2DistBlader + r3r2DistBlader);
    
    BladeCheckr = SegmentLengthBlader/BladeLength;
    
    count = 0;
    
    if (MinCheckRoot < RootCheckr) && (RootCheckr < MaxCheckRoot)
        count = count + 1;
    else
        fprintf('BAD /r/ Frame Root Displacement | Production: %i\n',ar)
    end
    
    if (MinCheckDorsum < DorsumCheckr) && ( DorsumCheckr < MaxCheckDorsum)
        count = count + 1;
    else
        fprintf('BAD /r/ Frame Dorsum Displacement | Production: %i\n',ar)
    end
    
    if (MinCheckBlade < BladeCheckr) && (BladeCheckr < MaxCheckBlade)
        count = count + 1;
    else
        fprintf('BAD /r/ Frame Blade Displacement | Production: %i\n',ar)
    end
    if count == 3
        fprintf('All /r/ Frame Displacements are GOOD\n')
    end
end