function compare_half_life(d1,d2);
% compare_half_life(d1,d2);
set(figure(5),'position',[561   544   422   405])
clf
plot( d1.halflife, d2.halflife, 'ko')
hold on
for m = 1:length( d1.halflife );
    plot( d1.halflife(m)*[1-d1.rel_error(m),1+d1.rel_error(m)],d2.halflife(m)*[1 1],'k-');
    plot( d1.halflife(m)*[1 1],d2.halflife(m)*[1-d2.rel_error(m),1+d2.rel_error(m)],'k-');
end

hold on
plot( [0 7], [0 7],'k');
xlim(([0 7 ])); ylim([0 7 ]);
xlabel( sprintf('half life (h), %s',d1.output_dir),'interp','none' ); 
ylabel( sprintf('half life (h), %s',d2.output_dir),'interp','none' ); 
set(gcf, 'PaperPositionMode','auto','color','white');