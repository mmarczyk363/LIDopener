function plot_enzymes(data_enz_norm,data_enz_canc,enzymes_expr,enzymes_expr_isos,sel_enzymes)

%isozymes plots for all enzymes
n = size(data_enz_norm,2);
colors = [0.8500,0.3250,0.0980;0,0.4470,0.7410;0.9290,0.6940,0.1250;0.4940,0.1840,0.5560;0.4660,0.6740,0.1880;0.3010,0.7450,0.9330;0.6350,0.0780,0.1840];
colors = repmat(colors,40,1);
p = zeros(size(sel_enzymes,1),1);
for a=1:size(sel_enzymes,1)
    ind1 = find(strcmp(enzymes_expr.ID,sel_enzymes.ID{a}));
    ind2 = find(enzymes_expr_isos.ind == ind1);
    
    %move target isozyme to the beginning
    genes_tmp = enzymes_expr_isos.gene(ind2);
    target_ind = strcmp(genes_tmp,sel_enzymes.GeneName{a});
    ind2 = [ind2(target_ind);ind2(~target_ind)];
    
    %get logFC
    data_test = data_enz_canc(ind2,:) - data_enz_norm(ind2,:); 
    if size(data_test,1) == 2
%         [~,p_canc] = ttest2(data_enz_canc(ind2(1),:),data_enz_canc(ind2(2),:));
        [~,p(a)] = ttest(data_test(1,:) - data_test(2,:));
    else
%         p_canc = anova1(data_enz_canc(ind2,:)',[],'off');
        p(a) = anova1(data_test',[],'off');
    end

    %expression difference between isozymes
    figure; hold on; box on
    boxplot(data_test','Labels',enzymes_expr_isos.gene(ind2)) 
    plot(xlim,[0,0],'k');
    tmp = axis(); text(tmp(1),tmp(3),['p = ' num2str(p(a))],'VerticalAlignment','bottom')
    title([enzymes_expr.ID{ind1},' - ',enzymes_expr.Name{ind1}])
    ylabel('Cancer - normal')
    saveas(gcf,['res/LID/enzymes_plots/enz_box_' enzymes_expr.ID{ind1} '.png']);
    saveas(gcf,['res/LID/enzymes_plots/enz_box_' enzymes_expr.ID{ind1} '.fig']);
    close gcf

    %expression diff. between conditions, separate isozymes, all
    figure; hold on; box on;
    for b=1:n
        data_plot = [data_enz_norm(ind2,b),data_enz_canc(ind2,b)];
        for c=1:size(data_plot,1)
            plot(data_plot(c,:),'Color',colors(c,:),'LineWidth',0.5,'Marker','.');
        end
    end
    err_n = 1.96*std(data_enz_norm(ind2,:)')/n;
    err_c = 1.96*std(data_enz_canc(ind2,:)')/n;
    errorbar([mean(data_enz_norm(ind2,:),2)';mean(data_enz_canc(ind2,:),2)'],[err_n;err_c],'LineWidth',2,'Color','k')
    legend(enzymes_expr_isos.gene(ind2),'Location','East')
    xlim([0.8,2.2])
    ylabel('Log2 expression');
    title([enzymes_expr.ID{ind1},' - ',enzymes_expr.Name{ind1}])
    tmp = axis(); text(tmp(1),tmp(3),['p = ' num2str(p(a))],'VerticalAlignment','bottom')
    set(gca,'XTick',1:2,'XTickLabel',{'Normal','Cancer'})
    % set transparency to lines
    lines = findobj(gcf, 'Type', 'Line');
    for l=1:length(lines); lines(l).Color(4) = 0.4; end
    saveas(gcf,['res/LID/enzymes_plots/enz_all_expr' enzymes_expr.ID{ind1} '.png']);
    saveas(gcf,['res/LID/enzymes_plots/enz_all_expr' enzymes_expr.ID{ind1} '.fig']);
    close gcf
end
