
clear all
close all
clc
Sign = dir(['...\' '/' '*.mat']);
NumSig = size(Sign,1);
ASF=[1:885];
ZZix=1;
for himland=1:5
    for i=1:23%NumSig
        feature1=[];
        Label=[];
        load(['...\',Sign(i).name])%Signal
        load(['...\',Sign(i).name])%Vigilance level
        zz=0;
        lim=1;
        final1=[];
        feature=[];
        for i=1:5
            feature=[feature;reshape(LOO(:,:,i),[17,885])];
        end
        feature1=[feature1,feature];
        Label=[Label,perclos'];
        X_t=feature1';
        clear feature1
        t=Label';
        y_t=t;
        X=X_t;
        X( (177*himland-176:177*himland),:)=[];
        feature1=X_t( (177*himland-176:177*himland),:);
        y=y_t;
        y( (177*himland-176:177*himland))=[];
        PERCLOS=y_t( (177*himland-176:177*himland));
        X_v=X;
        y_v=y;
        clear Label
        threshold=0.001;
        X_O=X;
        Y_O=feature1;
        iteration_pred4=[];
        for qiy=1:size(Y_O,1)
            X=X_v;
            y=y_v;
            for wiy=1:size(X,1)
                Dist(wiy)=norm(Y_O(qiy,:)-X (wiy,:));
            end
            [a,b]=find(Dist<5.*(((max(Dist)+min(Dist))/4)));
            X2=X(b,:);
            y2=y(b,:);
            X=X2;
            y=y2;
            clear b Dist X2 y2
            for i=1:size(X,1)
                for j=i+1:size(X,1)
                    K(i,j)=norm(X(i,:)-X(j,:));
                end
            end
            K=[K;zeros(1,size(K,2))];
            K=K+K';
            for i=1:size(feature1,1)
                for j=1:size(X,1)
                    kk1(i,j)=norm(Y_O(i,:)-X(j,:));
                end
            end
            X=K;
            feature1=kk1;
            clear K kk1
            M=1;
            for i=1:M
                w{i}=rand(size(X,2) ,1);
            end
            for i=1:M
                w{i}=w{i}./norm(w{i});
            end
            for i=1:M
                [g{i},~,~] = fit(X*w{i},y,'smoothingspline','SmoothingParam',0.01);
            end
            iteration_pred=0;
            for i=1:M
                iteration_pred=iteration_pred+g{i}(X*w{i});
            end
            init=iteration_pred;
            updated_error=norm(y-iteration_pred);
            uuo=1;
            wiz=w;
            giz=g;
            while (updated_error>threshold &uuo<500)
                for i=1:M
                    other_w=w;
                    other_w(i)=[];
                    other_g=g;
                    other_g(i)=[];
                    ppr_pred=0;
                    for j=1:size(other_w,2)
                        ppr_pred=ppr_pred+other_g{j}(X*other_w{j});
                    end
                    y_residual=y-ppr_pred;
                    w{i}=update_W(X,y_residual,w{i},g{i});
                    w{i}=w{i}./norm(w{i});
                    [g{i},~,~] = fit(X*w{i},y_residual,'smoothingspline','SmoothingParam',0.01);
                end
                iteration_pred=0;
                for i=1:M
                    iteration_pred=iteration_pred+g{i}(X*w{i});
                end
                updated_error=norm(y-iteration_pred);
                ll(uuo)=updated_error;
                uuo=uuo+1;
            end
            iteration_pred1=iteration_pred;
            iteration_pred3=0;
            for i=1:M
                iteration_pred3=iteration_pred3+g{i}(feature1(qiy,:)*w{i});
            end
            iteration_pred3(iteration_pred3>1)=1;
            iteration_pred3(iteration_pred3<0)=0;
            iteration_pred4=[iteration_pred4,iteration_pred3]
        end
        clf
        plot(iteration_pred4,'o')
        hold on
        plot(PERCLOS,'.')
        pause(0.2)
        RMSE(ZZix)=sqrt(sum((iteration_pred4'-PERCLOS).^2)./size(PERCLOS,1))
        ZZix=ZZix+1;
        clearvars -except ZZix RMSE Sign NumSig himland 
    end
end