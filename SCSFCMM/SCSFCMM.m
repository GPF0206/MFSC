function [clusterO,labelsO,U] = SCSFCMM(idx,idy,labdata,cluster_N,compactP,gridStepX,gridStepY,K,cluster,data,I)
invt=1/((gridStepX+gridStepY)/(compactP*2))^2;
GFSrange=8;
x11=zeros(idx*idy,1);
y11=zeros(idx*idy,1);
labelsO=zeros(idx*idy,1);
memb=zeros(idx*idy,GFSrange);
p=0;
q=2;
t=2;
U=zeros(idx*idy,GFSrange);
alex=0.1;
neighIndices = cell(idx*idy, 1);
for iter=1:10
    %初始化
    H=zeros(idx*idy,GFSrange);
    Lab_mat_c=zeros(idx*idy,GFSrange);
    countF=zeros(idx*idy,1);
    clusterS=zeros(size(cluster_N,2),1);
    %计算距离
    for i=1:size(cluster_N,2)
        %每个聚类中心的搜索区域
        x1=max(cluster_N(4,i)-gridStepX,1);
        y1=max(cluster_N(5,i)-gridStepX,1);
        x2=min(cluster_N(4,i)+gridStepY,idx);
        y2=min(cluster_N(5,i)+gridStepY,idy);
        for j=floor(x1):floor(x2)
            for k= floor(y1):floor(y2)
                index=(j-1)*idy+k;
                x11(index)=j+1;
                y11(index)=k+1;
                dist=(labdata(index,1)-cluster_N(1,i))^2+(labdata(index,2)-cluster_N(2,i))^2+(labdata(index,3)-cluster_N(3,i))^2;
                distxy=(j-cluster_N(4,i))^2+(k-cluster_N(5,i))^2;
                dist=dist+invt*distxy;
                dist=sqrt(dist);
                if countF(index)<GFSrange
                    countF(index)=countF(index)+1;
                    memb(index,countF(index))=i;
                    Lab_mat_c(index,countF(index))=dist;
                else
                    %找出距离最大的
                    maxDist=max(Lab_mat_c(index,:),2);
                    if dist<maxDist
                        loc=find(Lab_mat_c(index,:)==maxDist);
                        Lab_mat_c(index,loc)=dist;
                        memb(index,loc)=i;
                    end
                end
            end
        end
    end
    for i=1:idx*idy
        for j=1:countF(i)
            temp=0;
            for k=1:countF(i)
                temp=temp+(Lab_mat_c(i,j)/Lab_mat_c(i,k))^(2.0/(t-1));
            end
            U(i,j)=1/(temp+eps);
        end
    end

   for i = 1:idx*idy
    neighIndices{i} = Neighbourhood(idx, idy, i, 8);
   end
   for i = 1:idx*idy
    neigh = neighIndices{i};
    for j = 1:countF(i)
        temp = 0;
        for k = 1:numel(neigh)
            w = find(memb(neigh(k), :) == memb(i, j), 1);
            if ~isempty(w)
                temp = temp + U(neigh(k), w);
            end
        end
        H(i, j) = temp + U(i, j);
    end
   end
   U1=(U.^p).*(H.^q)./(sum((U.^p).*(H.^q),2)*ones(1,size(U,2)));
   U1_t_powers = U1.^2;
   
       %% CIELAB
    cluster_T=cluster_N(1:3,:)+eps;
    if iter~=1
        d_r=data(:,1);
        d_g=data(:,2);
        d_b=data(:,3);
        for i=1:max(labelsO(:))
            if sum(sum(labelsO==i))~=0
                l_r=median(d_r(labelsO==i));
                l_g=median(d_g(labelsO==i));
                l_b=median(d_b(labelsO==i));
                cluster_T(:,i)=[l_r l_g l_b]'+eps;
            end
        end
        cluster_T=colorspace('Lab<-RGB',cluster_T');
        cluster_T=cluster_T';
    end
    %%
   %更新C
    cluster_old=cluster_N;
    for i=1:idx*idy
        memb_i = memb(i, 1:countF(i));
        cluster_N(1:3, memb_i) = cluster_N(1:3, memb_i) + U1_t_powers(i, 1:countF(i)).* labdata(i, :)';
        cluster_N(4, memb_i) = cluster_N(4, memb_i) + U1_t_powers(i, 1:countF(i)).* x11(i);
        cluster_N(5, memb_i) = cluster_N(5, memb_i) + U1_t_powers(i, 1:countF(i)).* y11(i);
        for j=1:countF(i)
        clusterS(memb(i,j))=U1_t_powers(i,j)+clusterS(memb(i,j));
        end
        maxIndex=1;
        if countF(i)>1
            maxU=max(U(i,:));
            maxIndex=find(U(i,:)==maxU);
        end
        labelsO(i)=memb(i,maxIndex(1));
    end
   [~,Y,~,Cluster]= Self_SC1(cluster_T',K,cluster,3,1);
    %加入聚类中心
    domian0=zeros(5,K);
    domian1=sum(Y.^2,2);
    for i=1:cluster
        domian0(1:3,:)=domian0(1:3,:)+((Cluster(:,i)*ones(1,K)).*(ones(size(Cluster,1),1)*Y(:,i)').^2);
    end
    cluster_N=(cluster_N+alex*domian0)./((ones(size(cluster_N,1),1)*(clusterS)')+alex*(domian1.*ones(1,5))');
    if sum(sum(abs(cluster_old-cluster_N)))<=0.9||(iter==10)
        break;
    end
    for h = 1:K
            Label_seeds = DeserializeImg(labelsO,idx,idy);
            sxx=floor(cluster_N(4,h));
            syy=floor(cluster_N(5,h));
            if (sxx==0&&syy==0)
             continue;
            end
            % K-medoids重新更新超像素中心点
            [row1, col1] = find(Label_seeds == h);
            if ~isempty(row1) && ~isempty(col1)
                min_dist = Inf;
                min_row1 = 0;
                min_col1 = 0;
                for idx1 = 1:size(row1,1)
                    dc1 = (I(row1(idx1),col1(idx1),1)-cluster_N(1,h))^2+(I(row1(idx1),col1(idx1),2)-cluster_N(2,h))^2+(I(row1(idx1),col1(idx1),3)-cluster_N(3,h))^2;
                    ds1=(row1(idx1)-sxx)^2+(col1(idx1)-syy)^2;
                    d =sqrt(dc1+invt*ds1);
                    if d < min_dist
                        min_dist = d;
                        min_row1 = row1(idx1);
                        min_col1 = col1(idx1);
                    end
                end
                cluster_N(1,h) = I(min_row1, min_col1, 1);
                cluster_N(2,h) = I(min_row1, min_col1, 2);
                cluster_N(3,h) = I(min_row1, min_col1, 3);
            end
    end
end
    clusterO=cluster_N;
end

