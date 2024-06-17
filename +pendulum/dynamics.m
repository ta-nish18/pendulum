classdef dynamics < handle
    
    properties
        num         %ノードの個数
        m_vec       % n個のノードの質量のベクトル
        L_vec       % n個のノードの紐の長さのベクトル
        g     = 9.8 % 重力加速度


        % 参考元：https://note.com/sciencecafe_mc2/n/ne5077e5a62f1
        % 高速化のためθに関わらない行列部分のみ予め計算しておく
        Amat_temp
        Bmat_temp_diag
        Bmat_temp_nodiag
    end


    methods
        function obj = dynamics(num)
            obj.num = num;
        end
    
        function G = plot(obj,ax)
            if nargin < 2
                ax = gca;
            end

            
            s = 1:obj.num;
            e = 1+ s;
            G = plot(ax, graph(s,e) );

            G.MarkerSize = [1;obj.m_vec(:) *10];
            G.NodeLabel = ["", "m"+(1:obj.num)];
            G.EdgeLabel = "L"+(1:obj.num);
            G.NodeFontWeight = 'bold';
            G.EdgeFontWeight = 'bold';
            G.NodeFontSize = 20/sqrt(obj.num);
            G.EdgeFontSize = 20/sqrt(obj.num);


            c = colormap(ax);
            cidx = linspace(1,size(c,1),obj.num);
            G.NodeColor  = [ [0,0,0,]; c(round(cidx),:)];
            
            G.EdgeColor = 'k';

            hold(ax,'on')
            G.XData = G.XData - G.XData(1);
            G.YData = G.YData - G.YData(1);
            yline(ax, 0, 'k', 'LineWidth',4);
            hold(ax,'off')
        end

        % 微分方程式を構成するメソッド
            function out = dx(obj,x)
                theta     = x( 1:obj.num);
                velocity  = x((1:obj.num) + obj.num);
    
                diff_theta = theta - theta.';
    
                Amat = obj.Amat_temp .* cos(diff_theta);
                Bmat =  obj.Bmat_temp_diag   * diag( sin(theta) ) ...
                       +obj.Bmat_temp_nodiag * diag( velocity.^2) .* sin(diff_theta);
    
                dtheta    = velocity;
                dvelocity = - Amat \ Bmat * ones(obj.num,1);
    
                out = [dtheta;dvelocity];
            end
    
            function calc_ABmat(obj)
                if numel(obj.L_vec)~=obj.num || numel(obj.m_vec)~=obj.num
                    return
                end

                Mmat = zeros(obj.num);
                for i = 1:obj.num
                    mi = obj.m_vec(i);
                    Mmat(1:i,1:i) = Mmat(1:i,1:i) + mi * ones(i);
                end
    
                obj.Amat_temp = Mmat * diag(obj.L_vec);
                obj.Bmat_temp_diag   = diag(diag(Mmat)) * obj.g;
                obj.Bmat_temp_nodiag = (Mmat - diag(diag(Mmat)) ) * diag(obj.L_vec);
            end

        % データセットに関するメソッド
            function set.num(obj,num)
                obj.num = num;

                m = ones(1,num);
                n_min = min(numel(obj.m_vec),num); %#ok
                m(1:n_min) = obj.m_vec(1:n_min);   %#ok
                obj.m_vec = m;                     %#ok

                L = ones(1,num);
                n_min = min(numel(obj.L_vec),num); %#ok
                L(1:n_min) = obj.L_vec(1:n_min);   %#ok
                obj.L_vec = L;                     %#ok
            end
    
            function set.m_vec(obj,v)
                check(v,obj.num)%#ok
                obj.m_vec = v(:).';
                obj.calc_ABmat;
            end
    
            function set.L_vec(obj,v)
                check(v,obj.num)%#ok
                obj.L_vec = v(:).';
                obj.calc_ABmat;
            end
    end
end

function check(vec,num)
    if numel(vec) ~= num
        error('ベクトル要素数がnumと一致していません。')
    end
end