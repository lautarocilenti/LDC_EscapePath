function [] = PL_TestSearch(data)
%PL_TESTSEARCH 
PL_TestContour(data);
msLog = data.msLog;
sh = msLog{end}{3};
theta0 = msLog{1}{1};

myVideo = VideoWriter(sprintf("Videos/%s.mp4",data.M.searchAlgorithm),'MPEG-4'); %open video file
myVideo.FrameRate = 2;  %can adjust this, 5 - 10 works well for me
open(myVideo)
pauseTime = .5;
for i = 1:length(sh.thetaCurrent)
    x= -10; y = -10;
    plot(x,y,'ok','MarkerFaceColor','k','markersize',5)
    hold on
    quiver(x,y,x,y,'b','linestyle',':','linewidth',2);
    plot(x,y,'*r','MarkerFaceColor','r','markersize',12)
    legendstr = {"current step","step direction","new step"};
    legend(legendstr,'AutoUpdate','off')
    title(sprintf("%s",data.M.searchAlgorithm))
    
    PL_TestContour(data);
    axis([-5 5 -5 5])

    x = sh.thetaCurrent{i}(1,:);
    y = sh.thetaCurrent{i}(2,:);
    plot(x,y,'ok','MarkerFaceColor','k','markersize',5)
    pause(pauseTime)
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
  
    uI = unique(sh.iNewSearch{i});
    for j = 1:length(uI)
        dir = sh.stepDirectionSearch{i}(:,sh.iNewSearch{i} == uI(j));
        for k = 1:size(dir,2)
            if contains(data.M.searchAlgorithm,"Simplex")
                ii = 1:length(x);
                for p1 = ii
                    jj = find(ii~=p1);
                    xp = x(jj) -x(p1);
                    yp = y(jj) - y(p1);
                    for p2 = 1:length(jj)
                        quiver(x(p1),y(p1),xp(p2),yp(p2),'g','linestyle','--','linewidth',1.5,'ShowArrowHead',false);
                    end
                    
                end
                ii = find(ii~=uI(j));
                c1 = mean(x(ii));
                c2 = mean(y(ii));
                quiver(c1,c2,dir(1,k),dir(2,k),'xb','linestyle',':','linewidth',2);
            else
                quiver(x(uI(j)),y(uI(j)),dir(1,k),dir(2,k),'b','linestyle',':','linewidth',2);
            end
        end
        
    end
    frame = getframe(gcf); %get frame
        pause(pauseTime)
    writeVideo(myVideo, frame);
    xNew = sh.thetaNew{i}(1,:);
    yNew = sh.thetaNew{i}(2,:);
    plot(xNew,yNew,'*r','MarkerFaceColor','r','markersize',10,'linewidth',1.5)
    pause(pauseTime)
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    writeVideo(myVideo, frame);
    hold off
    
end

close(myVideo)

end

