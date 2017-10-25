function [mm_difference] = get_cochlear_distance(f1,f2)

%Chinch cochlear length
CochlearLength = 20.1;%mm. Based on Muller et al 2010. Hear Res. 268: 184-193.
d1 = absolutedistance(f1);%ipsi
d2 = absolutedistance(f2);%contra

mm_difference = d2-d1;%by this convention, ipsi-higher tuning (i.e., ipsi more basal) the difference in mm will be positive

    function [d] = absolutedistance(f)
        d = 61.2-42.2*log10(f);%Based on Muller et al 2010. Hear Res. 268: 184-193.
        d = (d/100)*CochlearLength;
    end
end