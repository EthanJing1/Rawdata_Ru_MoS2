function createaxes(Parent1, X1, YMatrix1)
%CREATEAXES(Parent1, X1, YMatrix1)
%  PARENT1:  axes parent
%  X1:  x 数据的向量
%  YMATRIX1:  y 数据的矩阵

%  由 MATLAB 于 14-Jun-2023 15:51:54 自动生成

% 创建 axes
axes1 = axes('Parent',Parent1,...
    'Position',[0.13 0.11 0.753928571428571 0.815]);
hold(axes1,'on');

% 使用 plot 的矩阵输入创建多行
plot1 = plot(X1,YMatrix1,'Parent',axes1);
set(plot1(4),'Color',[0 0 0]);
set(plot1(5),'Color',[0 0 0]);
set(plot1(6),'Color',[0 0 0]);

% 创建 xlabel
xlabel('Time');

Legend('1.71','1.74','1.78');

% 取消以下行的注释以保留坐标区的 Y 范围
% ylim(axes1,[-0.14319 0.011153]);
box(axes1,'on');
hold(axes1,'off');
% 设置其余坐标区属性
set(axes1,'ContextMenu');
