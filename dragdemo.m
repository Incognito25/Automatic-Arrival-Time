function dragdemo(x1,y1,x2,y2)
global x y uifig
% DRAGDEMO
%   This function is a demo for the draggable.m function and should be
%   distributed along with it. It basically presents some of draggable.m's
%   features: users of draggable.m are invited to read dragdemo's source.
%
%   'polymove'      A polygon with draggable vertices is drawn;
%                   draggable's "motionfcn" argument is used to redraw the
%                   polygon each time a vertex is moved.
%
% (C) Copyright 2004-2020
% François Bouffard
% fbouffard@gmail.com

% Modified by Silas Abraham 03-05-2022

% if nargin == 0
%     demotitle = 'polymove';
% end
% 
% switch lower(demotitle)
%     
%     case 'polymove'
        
        % POLYMOVE demo
        %   The 'motionfcn' argument of draggable.m is used here to
        %   redraw the polygon each time one of its vertex is moved.
        
        %   Setting up the figure and axes.
        
        uifig=figure;
        st=plot(x,y,'b');
        xlim([-5 50]);
        ylim([-1.1 1.1]);
        set(gcf, 'WindowState','maximized');

        hold off
%         set(gca,'Xlim',xlimit,'YLim',ylimit);
        
        % Creating the polygon vertices
        
        hold on;
        v1 = plot(x1(1),y1(1));
%         v2 = plot(x1(2),y1(2));
        v3 = plot(x2(1),y2(1));
%         v4 = plot(x2(2),y2(2));
        vv = [v1 v3];
        % For visible vertices:
        set(vv,'Marker','o','MarkerSize',7,'MarkerEdgeColor','r');
        
        % For invisible vertices:
        %set(vv,'Marker','o','MarkerSize',10,'MarkerFaceColor','none', ...
        %       'MarkerEdgeColor','none');
        
        % Saving the vertex vector as application data
        % in the current axes (along with empty element p which will
        % later hold the handle to the polygon itself)
        setappdata(gca,'st',st);
        setappdata(gca,'vv',vv);
        setappdata(gca,'p',[]);
%         setappdata(gca,'a',[]);
        % Calling draggable on each of the vertices, passing as an
        % argument the handle to the redraw_poly fucntion (see below)
        
        draggable(v1,@redraw_poly)
%         draggable(v2,@redraw_poly);
        draggable(v3,@redraw_poly);
%         draggable(v4,@redraw_poly);
        
        % Finally we draw the polygon itself using the redraw_poly
        % function, which can be found below
        
        redraw_poly;
        
% end


% -----------------------------------------------------------------------
% Function REDRAW_POLY
%   This function is passed as the 'motionfcn' argument to draggable.m in
%   the POLYMOVE demo. It recieves the handle to the object being dragged
%   as its only argument, but it is not actually used in this function.

function redraw_poly(h)
global xdata ydata
% Deleting the previous polygon

delete(getappdata(gca,'p'));
% delete(getappdata(gca,'a'));

% Retrieving the vertex vector and corresponding xdata and ydata
st = getappdata(gca,'st');
vv = getappdata(gca,'vv');
xdata = cell2mat(get(vv,'xdata'));
ydata = cell2mat(get(vv,'ydata'));
% Plotting the new polygon and saving its handle as application data

% p1 = plot([xdata(1) xdata(2)],[ydata(1) ydata(2)],'--r');
% p2 = plot([xdata(3) xdata(4)],[ydata(3) ydata(4)],'--r');

p=rectangle('Position',[xdata(1),ydata(1),abs(xdata(2)-xdata(1)),abs(ydata(2)-ydata(1))]);

% p=[p1 p2];
setappdata(gca,'p',p);

% Putting the vertices on top of the polygon so that they are easier
% to drag (or else, the polygone line get in the way)

set(gca,'Children',[vv p st]);

