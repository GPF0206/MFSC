clear all
close all
clc

%% 输入图片
%% 超像素输入
I0=imread('108_0878.jpg'); 
I=I0;
[idx,idy,~]=size(I);
I=rgb2lab(I);
length=idx*idy;
N=length;%
K=256;%k个种子点，聚类中心
m=25;%调节系数
S=floor(sqrt(N/K));%超像素之间的距离Mend
grid_stepX=round(size(I,1)/sqrt(K));
grid_stepY=round(size(I,2)/sqrt(K));
% [Labels, numlabels,x,y,l,a,b] = FuzzySLIC(I0,K,25,1,0.2);
% seeds = [l';a';b';y';x'];
[~,seeds]=SlicF(I,K,grid_stepX,grid_stepY );
%% 绘图
cluster=3;
g=3.95;
[seeds,labelsByFuzzySLIC,U]=SCSFCMM(idx,idy,SerializeImg(I),seeds',25,grid_stepX,grid_stepY,K,cluster,SerializeImg(I0),rgb2lab(I0));
Labels=DeserializeImg(labelsByFuzzySLIC,idx,idy);
%% 增强连通性
S_Search=16;%合并孤立点的筛选面积
%将小区域合并
% LabelsOld = Labels;
for i = 1:K
    %寻找第i标签superpixels标注
    emptylabels = zeros(idx,idy);
    %将每一类的所有点标注
    emptylabels(Labels == i) = 1;
    %找到所有隶属于这一超像素的像素，然后置1
    [L_tmp,num] = bwlabel(emptylabels,4);
    for j=1:num
        %将某超像素里面的像素个数小于预定个数的超像素进行拆分处理
        if sum(sum(L_tmp==j))<S_Search
            [x,y]=find(L_tmp==j);
            A = [];
            for ix=1:1:size(x)
                up=Labels(x(ix),min(y(ix)+1,idy));
                down=Labels(x(ix),max(y(ix)-1,1));
                right=Labels(min(x(ix)+1,idx),y(ix));
                left=Labels(max(x(ix)-1,1),y(ix));
                if up~=Labels(x(ix),y(ix))
                    A = [A, up];
                end
                if down~=Labels(x(ix),y(ix))
                    A = [A, down];
                end
                if left~=Labels(x(ix),y(ix))
                    A = [A, left];
                end
                if right~=Labels(x(ix),y(ix))
                    A = [A, right];
                end
            end
            if size(A)~=0
                Labels(x,y) = mode(A,2);
            end
        end
    end
end
%% 绘图
figure
imshow(lab2rgb(I),[]);hold on;
plot(seeds(:,1),seeds(:,2),'.r');hold on;
contour(Labels,K,'black');
[~,center_p,Num_p,center_lab]=Label_image(I0,Labels);
%% Self_SC聚类
X = [center_lab,center_p];
[Label,~,~,~]= Self_SC1(X,K,cluster,3,0.69);
Lr2=zeros(size(Labels,1),size(Labels,2));
for iter=1:max(Labels(:))
    Lr2=Lr2+(Labels==iter)*Label(iter);
end
Lseg=Label_image(I0,Lr2);
figure,imshow(Lseg);