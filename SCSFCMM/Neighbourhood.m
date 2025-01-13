function [label,distance]=Neighbourhood(idx,idy,num,nums)
if num<1||num>idx*idy
    error('input num error');
end
if nums~=4&&nums~=8
    error('input nums error');
end

label=[];
distance=[];
if nums==4
    if mod(num,idy)~=0
        label=[label,num+1];
        distance=[distance,1];
    end
    if mod(num,idy)~=1
        label=[label,num-1];
        distance=[distance,1];
    end
    if ceil(num/idy)~=idx
        label=[label,num+idy];
        distance=[distance,1];
    end
    if ceil(num/idy)~=1
        label=[label,num-idy];
        distance=[distance,1];
    end
    return
end
if nums==8
    if mod(num,idy)~=0
        label=[label,num+1];
        distance=[distance,1];
    end
    if mod(num,idy)~=1
        label=[label,num-1];
        distance=[distance,1];
    end
    if ceil(num/idy)~=idx
        label=[label,num+idy];
        distance=[distance,1];
    end
    if ceil(num/idy)~=1
        label=[label,num-idy];
        distance=[distance,1];
    end
    if mod(num,idy)~=0&&ceil(num/idy)~=idx
        label=[label,num+1+idy];
        distance=[distance,sqrt(2)];
    end
    if mod(num,idy)~=0&&ceil(num/idy)~=1
        label=[label,num+1-idy];
        distance=[distance,sqrt(2)];
    end
    if mod(num,idy)~=1&&ceil(num/idy)~=idx
        label=[label,num-1+idy];
        distance=[distance,sqrt(2)];
    end
    if mod(num,idy)~=1&&ceil(num/idy)~=1
        label=[label,num-1-idy];
        distance=[distance,sqrt(2)];
    end
    return
end