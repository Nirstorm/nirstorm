function segmentation_map = nst_prepare_segmentation(segmentation_map,varagin)
%NST_PREPARE_SEGMENTATION Prepare segmentation volumes for usages in
%NIRSTORM. Take in input as segmentation map and the value assigned to the
% different tissue and output the same map with the following label : 
% 1: skin,  2: skull, 3: CSF, 4: GM, 5: WM
skin=varagin{1};
skull=varagin{2};
CSF=varagin{3};
GM=varagin{4};
WM=varagin{5};

    if ~(skin==1 && skull==2  && CSF==3 && GM==4 && WM==5)

        old_map=segmentation_map;
        segmentation_map=zeros(size(segmentation_map));
        segmentation_map(old_map==skin)=1;
        segmentation_map(old_map==skull)=2;
        segmentation_map(old_map==CSF)=3;
        segmentation_map(old_map==GM)=4;
        segmentation_map(old_map==WM)=5;
    end

end

