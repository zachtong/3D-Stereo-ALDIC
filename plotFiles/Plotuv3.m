function Plotuv3(U,x,y)
%FUNCTION Plotuv(U,x,y )
% To plot DIC solved displacement fields
% ----------------------------------------------
%
%   INPUT: U                 Displacement cell 
%                            
%          x,y               x- and y-coordinates
%
%   OUTPUT: Plots of displacement fields.
%
% ----------------------------------------------
% Author: Jin Yang.  
% Contact and support: jyang526@wisc.edu -or- aldicdvc@gmail.com
% Last time updated: 11/2020.
% ==============================================


figure; surf(x,y,u);  set(gca,'fontSize',18); 
title('$x-$displacement $u$','FontWeight','Normal','Interpreter','latex');
axis tight; % set(gca,'XTick',[] );
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
box on; colormap jet;
 

figure; surf(x,y,v);  set(gca,'fontSize',18);  
title('$y-$displacement $v$','FontWeight','Normal','Interpreter','latex');
axis tight; % set(gca,'XTick',[] );
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
box on; colormap jet;

figure; surf(x,y,  w);  set(gca,'fontSize',18);  
title('$z-$displacement $w$','FontWeight','Normal','Interpreter','latex');
axis tight; % set(gca,'XTick',[] );
xlabel('$x$ (pixels)','Interpreter','latex'); ylabel('$y$ (pixels)','Interpreter','latex');
set(gcf,'color','w');
a = gca; a.TickLabelInterpreter = 'latex';
b = colorbar; b.TickLabelInterpreter = 'latex';
box on; colormap jet;

%% Some old codes
% print -painters -dpng -r600  fig_Sample14_Local_u.png
 
% figure; mesh(x,y,reshape( Exact(coordinatesFEM(:,1),2),M,N)); set(gca,'fontSize',20); view(-5, 30);title('Exact x displacement');
% axis(l);
% % print -painters -dpng -r600  fig_Sample14_Local_uExact.png
% 
% 
% figure; mesh(x,y, v);  set(gca,'fontSize',20); view(-5, 30);title('Computed y displacement'); axis tight
% l = axis;
% % print -painters -dpng -r600  fig_Sample14_Local_v.png
% 
% 
% figure; mesh(x,y,reshape(0* Exact(coordinatesFEM(:,1),2),M,N)); set(gca,'fontSize',20); view(-5, 30);title('Exact y displacement');
% axis(l);
% % print -painters -dpng -r600  fig_Sample14_Local_vExact.png


