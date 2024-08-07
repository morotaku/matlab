Red=[1 1 1];
for i=0:50
    Red=[Red;[1 1-0.02*i 1-0.02*i]];
end
for i=0:10
    Red=[Red;1 0 0];
end

%mycolors = [1 1 1; 1 . 0];
colormap(Red)