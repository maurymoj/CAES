function applystyle2plot(varargin)
% PURPOSE:  Applies figure style to the corrent figure. It changes the 
% properties of the current plot such as
% background color, font, font size, figure aspecct ratio, etc,
% automatically (no need to change those elements manually)
%---------------------------------------------------
% USAGE: applystyle2plots()
%        applystyle2plots('property',property_value)
%---------------------------------------------------
% RETURNS:
% --------------------------------------------------
% EXAMPLE:
% x=linspace(0,4*pi,100);
% y=cos(x);
% z=sin(x);
% plot(x,y,'ro-',x,z,'g--');
% title('The title goes here');
% xlabel('The x label goes here');
% ylabel('The y label goes here');
% legend('cos(x)','sin(x)')
% applystyle2plot()
% --------------------------------------------------
% SEE ALSO:
%---------------------------------------------------
% REFERENCES: 
%---------------------------------------------------
% REMARKS: May not work for more complicated types of plots
%---------------------------------------------------
% WRITTEN BY: Hugo T. C. Pedro, 2011 
%             Maury M. Oliveira Jr, 2018

    %sets the predefined styles
    Font = 'Arial';
    LabelSize = 14;
    AxisFontSize = 12;
    LegendSize = 12;
    AspectRatio = 1;
    MarkerSize = 4;
    LineWidth = 1;
    Margins = [0 0 0 0];
    PrintSize = 5;
    colormap('gray');
    
    %parses input arguments
    if nargin > 0    
        for i=1:2:nargin
            arg1 =  varargin{i};
            arg2 =  varargin{i+1};
            switch arg1
                case 'Font'
                    Font = arg2;
                case 'LabelSize'
                    LabelSize = arg2;
                case 'AspectRatio'
                    AspectRatio = arg2;
                case 'Margins'
                    Margins = arg2;
                case 'MarkerSize'
                    MarkerSize = arg2;
                case 'LegendSize'
                    LegendSize = arg2;     
                case 'AxisFontSize'
                    AxisFontSize = arg2;   
                case 'PrintSize'
                    PrintSize = arg2; 
                case 'LineWidth'
                    LineWidth = arg2; 
                otherwise
                    disp('Unknown argument');
                    return;
            end
        end 
    end

    % gets the handle of the current figure    
    hParent = gcf;
    
    %gets handles from children objects
    %lines
    hLine = findobj(hParent,'Type','line');
    %axes (not legends)
    hAxes = findobj(hParent,'Type','axes','-not', 'Tag','legend');
    %legends
    hLegend = findobj(hParent,'Tag','legend');
    %annotations
    hAnnotation = findobj(hParent,'Type','text');
    
    
    %changes the properties
    %sets the background color to white
    set(hParent,'Color',[1 1 1]);
        
    %sets the line style    
    set(hLine,'LineWidth',LineWidth,'MarkerSize',MarkerSize);
    %sets the marker color equal to the line color    
    lineColor = get(hLine,'Color');
    if (iscell(hLine))
        hLine = cell2mat(hLine);
    end    
    if (iscell(lineColor))
        lineColor = cell2mat(lineColor);
    end    
    for i=1:length(hLine)
        set(hLine(i),'MarkerFaceColor',lineColor(i,:)); 
    end

    %sets the axes style    
    set(hAxes,'FontName',Font,'FontSize',AxisFontSize,'Box','on','LineWidth',1);
    
    %sets the legend style    
    hLegendText = findobj(hLegend,'Type','text');
    set(hLegendText,'FontName',Font,'FontSize',LegendSize);

    %sets the text anotations style    
    set(hAnnotation,'FontName',Font,'FontSize',LegendSize);

    %sets the labels
    hXlabel = get(hAxes,'XLabel');
    if (iscell(hXlabel))
        hXlabel = cell2mat(hXlabel);
    end
    for i=1:length(hXlabel)
        set(hXlabel(i),'FontName',Font,'FontSize',LabelSize); 
    end

    hYlabel = get(hAxes,'YLabel');
    if (iscell(hYlabel))
        hYlabel = cell2mat(hYlabel);
    end
    for i=1:length(hYlabel)
        set(hYlabel(i),'FontName',Font,'FontSize',LabelSize);
    end
    
    %sets the title style    
    hTitle = get(hAxes,'Title');
    if (iscell(hTitle))
        hTitle = cell2mat(hTitle);
    end
    for i=1:length(hTitle)
        set(hTitle(i),'FontName',Font,'FontSize',LabelSize); 
    end
    
 
    %sets the image size and Aspect ratio
    if AspectRatio>=1
        pSize = [AspectRatio*PrintSize PrintSize];
    else
        pSize = [PrintSize PrintSize/AspectRatio];
    end

    set(hParent, 'PaperSize',pSize);
    set(hParent, 'PaperUnits', 'inches');
    set(hParent, 'PaperPosition', [Margins(1) Margins(2) pSize(1)-(Margins(1)+Margins(3)) pSize(2)-(Margins(2)+Margins(4))]);    
     
end
    