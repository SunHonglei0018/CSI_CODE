function arr = ThreeToOne(par)
temp(:,:) = [par(1,1,:);par(1,2,:);par(1,3,:)];
arr(:,1) = [temp(1,:),temp(2,:),temp(3,:)];
arr=arr';
end
