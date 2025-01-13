function [ out_img ] = DeserializeImg( input_img,idx,idy )
if size(input_img,2)==1
    out_img=zeros(idx,idy);
    for i=1:idx
        out_img(i,:)=input_img((i-1)*idy+1:i*idy);
    end
else
    out_img=zeros(idx,idy,size(input_img,2));
    for i=1:idx
        out_img(i,:,:)=input_img((i-1)*idy+1:i*idy,:);
    end
end
