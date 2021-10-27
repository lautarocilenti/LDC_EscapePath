function [] = PL_TestSearch(data)
%PL_TESTSEARCH 
PL_TestContour(data);
msLog = data.msLog;
sh = msLog{end}{3};
theta0 = msLog{1}{1};

myVideo = VideoWriter(sprintf("Videos/%s.avi",data.M.searchAlgorithm)); %open video file
myVideo.FrameRate = 3;  %can adjust this, 5 - 10 works well for me
open(myVideo)
pauseTime = .1;
for i = 1:length(sh.thetaCurrent)
    x= -10; y = -10;
    plot(x,y,'ok','MarkerFaceColor','k','markersize',3)
    hold on
    quiver(x,y,x,y,'b','linestyle',':','linewidth',1);
    plot(x,y,'*r','MarkerFaceColor','r','markersize',8)
    legendstr = {"current step","step direction","new step"};
    legend(legendstr,'AutoUpdate','off')
    title(sprintf("%s",data.M.searchAlgorithm))
    
    PL_TestContour(data);
    axis([-5 5 -5 5])

    x = sh.thetaCurrent{i}(1,:);
    y = sh.thetaCurrent{i}(2,:);
    plot(x,y,'ok','MarkerFaceColor','k','markersize',3)
    pause(pauseTime)
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
  
    uI = unique(sh.iNewSearch{i});
    for j = 1:length(uI)
        dir = sh.stepDirectionSearch{i}(:,sh.iNewSearch{i} == uI(j));
        for k = 1:size(dir,2)
            if contains(data.M.searchAlgorithm,"Simplex")
                ii = 1:3;
                ii = find(ii~=uI(j));
                c1 = mean(x(ii));
                c2 = mean(y(ii));
                quiver(c1,c2,dir(1,k),dir(2,k),'xb','linestyle',':','linewidth',1);
            else
                quiver(x(uI(j)),y(uI(j)),dir(1,k),dir(2,k),'b','linestyle',':','linewidth',1);
            end
        end
        
    end
    frame = getframe(gcf); %get frame
        pause(pauseTime)
    writeVideo(myVideo, frame);
    xNew = sh.thetaNew{i}(1,:);
    yNew = sh.thetaNew{i}(2,:);
    plot(xNew,yNew,'*r','MarkerFaceColor','r','markersize',8)
    pause(pauseTime)
    frame = getframe(gcf); %get frame
    writeVideo(myVideo, frame);
    writeVideo(myVideo, frame);
    hold off
    
end

close(myVideo)

end

