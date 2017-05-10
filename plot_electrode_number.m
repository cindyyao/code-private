function plot_electrode_number(position, rotation)

for i = 1:length(position)
    [t1 t2]=cart2pol(position(i,1),position(i,2));
    t1=t1+(rotation/(180/pi));
    [position(i,1),position(i,2)]=pol2cart(t1, t2);
end

for ee=1:size(position,1)
    labeltext = num2str(ee);
    text(position(ee,1),position(ee,2),labeltext,'Color', 'k','FontSize',10,...
        'HorizontalAlignment','Center','VerticalAlignment','Bottom');
    hold on
end
xlim([-400 400])
ylim([-400 400])
axis equal

end