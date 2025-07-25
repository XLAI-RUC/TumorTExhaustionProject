function dy = ObjFun(data,soly)
dy=sum(((data(:,2)-soly)/max(data(:,2))).^2);
end