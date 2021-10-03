function [halflife_mean,halflife_err] = half_life_compilation_plot(sample_names, sample_colors, outputs)
% [halflife_mean,halflife_err] = half_life_compilation_plot(sample_names, sample_colors, outputs)

markers = {'s','o','v','^','x','*','d','p'};

set(gcf,'pos',[64   579   800   376]);
clf
halflife_mean = nan * ones(1,length(sample_names));
halflife_err  = nan * ones(1,length(sample_names));
for i = 1:length(sample_names)
    sample_name = sample_names{i};
    h = []; h_err = [];
    for n = 1:length(outputs)
        output = outputs{n};
        idx = i;
        if ( idx > length(output.halflife)) idx = []; end;
%         if isfield( output, 'sample_names' )
%             idx = find(strcmp(output.sample_names,sample_name))
%         end
        if ~isnan( output.halflife(i) )
            h = [h, output.halflife(i)];
            h_err = [h_err, output.halflife(i) * output.rel_error(i)];
        end
    end
    halflife_mean(i) = sum( h./h_err.^2 )/sum(1./h_err.^2);

    % propagate bootstrap errors
    %halflife_err(i)  = sqrt( 1./sum(1./h_err.^2) );

    % std err = std deviation /sqrt(N)
    halflife_err(i)  = std(h)/sqrt(length(h));
end


% to set up for legends
for n = 1:length(outputs)
    plot( 1, 0.5, 'ko','markerfacecolor','k','markersize',5,'marker',markers{n} ); hold on;
    legend_titles{n} = outputs{n}.output_dir;
end


halflife_ref = halflife_mean(1);
for i = 1:length(sample_names)
    z = barh( i,halflife_mean(i)/halflife_ref,'facecolor',sample_colors{i} ,'edgecolor','none');
    hold on
    plot( halflife_mean(i)/halflife_ref + [-halflife_err(i) halflife_err(i)]/halflife_ref, i*[1 1], 'k','linew',2 ); hold on;
end

for n = 1:length(outputs)
    output = outputs{n};
    halflife = output.halflife;
    rel_error = output.rel_error;
    for i = 1:length(halflife)
        plot( halflife(i)/halflife_ref, i+0.1*randn(1), 'k','markerfacecolor','k','markersize',5,'marker',markers{n} ); hold on;
    end
end

plot( [ 1 1],0.5+[-1 length(sample_names)+1],'color',[0.5 0.5 0.5],'linew',2); hold on
ylim(0.5+[0 length(sample_names)])
set(gca,'ytick',[1:length(sample_names)],...
    'yticklabel', sample_names,...
    'ticklabelinterp','none',...
    'ydir','reverse');
set(gca,'fontweight','bold');
xlabel( 'Half life (rel. to reference)' );
set(gcf, 'PaperPositionMode','auto','color','white');
h=legend( legend_titles );
set(h,'interp','none','location','northoutside');
