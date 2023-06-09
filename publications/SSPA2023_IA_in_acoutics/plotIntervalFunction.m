function [f_smp,f_int] = plotIntervalFunction(f,a_int,b_int,varargin)

    % Calculate function value
    f_int = f(a_int,b_int);

    % Plot
    hold on
    if ~isempty(varargin)
        fsurf(f,varargin{1},'MeshDensity',10,'DisplayName','function')
        a_lim = varargin{1}(1:2);
        b_lim = varargin{1}(3:4);
    else
        fsurf(f,'MeshDensity',10,'DisplayName','function') 
        a_lim = xlim();
        b_lim = ylim();
    end
    f_lim = zlim();
    plot3(a_int.Bounds,[b_lim(1) b_lim(1)],[f_lim(1) f_lim(1)],...
                                                'r-','LineWidth',4,...
                                        'HandleVisibility','off')
    plot3([a_lim(1) a_lim(1)],b_int.Bounds,[f_lim(1) f_lim(1)],...
                                                'r-','LineWidth',4,...
                                        'HandleVisibility','off')
    h_f = fsurf(f,[a_int.Bounds b_int.Bounds],'MeshDensity',10,...
                                                'LineWidth',4,...
                                                'EdgeColor','r',...
                                        'HandleVisibility','off');
    f_smp = [min(h_f.ZData),max(h_f.ZData)];
    plot3([a_lim(1) a_lim(1)],[b_lim(2) b_lim(2)],f_smp,...
                                                'r-','LineWidth',4,...
                                        'DisplayName','valid')
    plot3([a_lim(1) a_lim(1)],[b_lim(2) b_lim(2)],f_int.Bounds,...
                                                'r:','LineWidth',4,...
                                        'DisplayName','interval')
end