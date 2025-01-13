function out_img = SerializeImg(img)
%%out_img类型为double
idx=size(img,1);
idy=size(img,2);
%声明为zeros为double类型
out_img=zeros(idx*idy,size(img,3));
if size(img,3)==1
    %gray image
    for i=1:idx
        out_img((i-1)*idy+1:(i-1)*idy+idy)=img(i,:);
    end
else
    %color
    for i=1:idx
        out_img((i-1)*idy+1:(i-1)*idy+idy,:)=squeeze(img(i,:,:));
    end
end
end

