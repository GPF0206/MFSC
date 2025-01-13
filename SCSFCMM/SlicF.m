function   [Labels,seeds]=SlicF(I,K,S1,S2 )
%I=imread('12003.jpg');
[idx,idy,c]=size(I);
% I=rgb2lab(I);
length=idx*idy;
N=length;%
%K=256;%k个种子点，聚类中心
m=25;%调节系数
S=floor(sqrt(N/K));%超像素之间的距离Mend
%超像素之间的横向距离
% S1=floor(idx/sqrt(K));
% %超像素之间的纵向距离
% S2=floor(idy/sqrt(K));

NewCluster=zeros(K,4);
% seeds=zeros(K,3);

%% 初始化种子点
% SeedVector=zeros(K,2);%存放种子点的坐标
XOffset=0;
YOffset=0;
n=floor(sqrt(K));
SeedVectorOld = zeros(K, 5);%labxy
% 种子点都为初始化超像素的中心
for j=1:1:n
    for i=1:1:n
        SeedVectorOld(i+n*(j-1),4)=min(floor(S1/2)+XOffset,idx);
        SeedVectorOld(i+n*(j-1),5)=min(floor(S2/2)+YOffset,idy);
        YOffset=YOffset+S2;
    end
    XOffset=XOffset+S1;
    YOffset=0;
end 
%% 3*3邻域选择最小梯度
G=I;
SeedVector=SeedVectorOld;

for i=1:1:idx-1
   for j=1:1:idy-1
        dx=I(i+1,j)-I(i,j);
        dy=I(i,j+1)-I(i,j);
        G(i,j)=dx+dy;
   end
end
%遍历所有种子点，移动到邻域内的最小的点的位置
%聚类中心 
for i=1:1:K
   Gx=SeedVector(i,4);
   Gy=SeedVector(i,5);
   A=G(max(Gx-1,1):min(Gx+1,idx), max(Gy-1,1):min(Gy+1,idy));
   [Xindex,Yindex]=find(A==min(min(A)));
   SeedVector(i,4)=SeedVector(i,4)+Xindex(1)-2;
   SeedVector(i,5)=SeedVector(i,5)+Yindex(1)-2;
end
seeds=SeedVector;
%% 初始化迭代变量
Labels = -1*ones(idx,idy);
Distance = inf(idx,idy);
error=0.0001;
residual=1;
I=double(I);
iternum = 1;
iters=1000;
%%
% while residual>error
while iternum<iters
 %% 对所有像素按照种子距离打上标记（进行分类）
 % 隶属度
    for k=1:1:K
        sx=floor(seeds(k,4));
        sy=floor(seeds(k,5));
        if (sx==0&&sy==0)
            continue;
        end
        % 2S*2S
        for i=max(1, floor(sx-floor(S1))+1):1:min(idx, floor(sx+floor(S1)))
            for j=max(1, floor(sy-floor(S2))+1):1:min(idy, floor(sy+floor(S2)))
                if sx >= 1 && sx <= idx && sy >= 1 && sy <= idy
                 dc = sqrt((I(i,j,1)-I(sx,sy,1))^2+(I(i,j,2)-I(sx,sy,2))^2+(I(i,j,3)-I(sx,sy,3))^2);
                 ds=sqrt((i-sx)^2+(j-sy)^2);
                end
%                 dc=sqrt((I(i,j,1)-I(sx,sy,1))^2+(I(i,j,2)-I(sx,sy,2))^2+(I(i,j,3)-I(sx,sy,3))^2);
%                 %当前像素和seed的欧氏距离
%                 ds=sqrt((i-sx)^2+(j-sy)^2);
                D= sqrt(dc*dc+m^2*(ds/S)^2);
%                 D= sqrt(dc*dc+m^2*(ds*0.1)^2);
                if D<Distance(i,j)
                   Distance(i,j) = D;
                   Labels(i,j) = k;
                end
            end
        end
    end
%% 迭代聚类中心
    %计算新的聚类中心
    NewCluster=zeros(K,6);
     for ix = 1:idx
            for iy = 1:idy
                label = Labels(ix,iy);
                if(label==-1)
                    continue
                end
%                 disp([label, ix, iy])
                %将当前类别号在x,y点的I,x,y值进行累加至NewCluster
%                 NewCluster(label,1) = NewCluster(label,1)+ix;%一类的x坐标之和
%                 NewCluster(label,2) = NewCluster(label,2)+iy;%一类的y坐标之和
%                 NewCluster(label,3) = NewCluster(label,3)+I(ix,iy,1)+I(ix,iy,2)+I(ix,iy,3);%灰度值之和
%                 NewCluster(label,4) = NewCluster(label,4)+1;%一类的总个数
                NewCluster(label,1) = NewCluster(label,1)+I(ix,iy,1);%l
                NewCluster(label,2) = NewCluster(label,2)+I(ix,iy,2);%a
                NewCluster(label,3) = NewCluster(label,3)+I(ix,iy,3);%b
                NewCluster(label,4) = NewCluster(label,4)+ix;%一类的x坐标之和
                NewCluster(label,5) = NewCluster(label,5)+iy;%一类的y坐标之和
                NewCluster(label,6) = NewCluster(label,6)+1;%一类的总个数
            end
     end
        SeedsOld = seeds;
        %求平均，重新计算该种子点的中心点
 %% 迭代种子点
        for i = 1:K
            seeds(i,1:5) = round(NewCluster(i,1:5)/(NewCluster(i,6)+eps));
        end
 %%
%         %判断是否收敛
        curErr = norm(norm(SeedsOld-seeds));
        if curErr<error
            break;
        end
    iternum=iternum+1;
end
disp(iternum)
%% 增强连通性
S_Search=16;%合并孤立点的筛选面积
%将小区域合并
LabelsOld = Labels;
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
% figure
% imshow(lab2rgb(I),[]);hold on;
% % plot(seeds(:,1),seeds(:,2),'.r');hold on;
% contour(Labels,K,'black');

end

