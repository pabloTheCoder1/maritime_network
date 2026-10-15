clear; clc; close all;

%costruzione della rete
% legge la lista di archi pesati dal file edges.csv
% (porti come nodi, rotte commerciali come archi, peso = intensità traffico)

T = readtable('edges.csv');
G = graph(T.source, T.target, T.weight);
N = numnodes(G);

% Visualizzazione della rete con layout a forza
figure('Name','Rete del trasporto marittimo globale');
p = plot(G, 'Layout','force', 'WeightEffect','direct');
p.LineWidth = 0.5 + 2.5 * G.Edges.Weight / max(G.Edges.Weight);
p.MarkerSize = 5;
title('Rete del trasporto marittimo globale');

% Evidenzia i nodi collo di bottiglia in rosso
bottlenecks = {'Suez','Panama','Malacca'};
highlight(p, bottlenecks, 'NodeColor','r', 'MarkerSize', 9);

% Grado, average path length, diametro, clustering coefficient
deg = degree(G);

D_unw = distances(G, 'Method','unweighted');
%matrice 31×31 dove la cella (i,j) contiene il numero di hop tra il nodo i e il nodo j
mask  = ~isinf(D_unw) & D_unw > 0;  %matrice booleana che pone a 0 celle infinite e <0
avg_path = mean(D_unw(mask));
diameter = max(D_unw(mask));
C = local_clustering_coefficient(G);

% simulazione di robustezza
% confronto tra rimozione mirata (hub) e rimozione casuale
% metriche misurate ad ogni step
% S(f): dimensione relativa della Giant Component
% l(f): APL residua sulla Giant Component

n_trials = 30;        % n trial per la rimozione casuale
steps = 0 : floor(N * 0.55);     % fino al 55% dei nodi rimossi
frac = steps / N;

% Rimozione mirata 
[~, ord] = sort(deg, 'descend'); %indici ordinati
hub_order = G.Nodes.Name(ord); %nodi ordinati per grado

gc_targ  = zeros(size(steps));
apl_targ = zeros(size(steps));
Gx = G;

for s = 1 : length(steps)
    if s > 1
        Gx = rmnode(Gx, hub_order{s-1});
    end
    [gc_targ(s), apl_targ(s)] = network_metrics(Gx, N);
end

% --- Rimozione casuale (media su n_trials) ---
rng(42);
gc_rand_acc  = zeros(size(steps));
apl_rand_acc = zeros(size(steps));

for t = 1 : n_trials
    perm = G.Nodes.Name(randperm(N));
    Gx   = G;
    for s = 1 : length(steps)
        if s > 1
            Gx = rmnode(Gx, perm{s-1});
        end
        [gc, apl] = network_metrics(Gx, N);
        gc_rand_acc(s)  = gc_rand_acc(s)  + gc;
        apl_rand_acc(s) = apl_rand_acc(s) + apl;
    end
end

gc_rand  = gc_rand_acc  / n_trials;
apl_rand = apl_rand_acc / n_trials;

% Plot robustezza 
figure('Name','Simulazione di robustezza');

subplot(1,2,1);
plot(frac, gc_targ, 'r-o', 'LineWidth',1.8, 'MarkerSize',4); hold on;
plot(frac, gc_rand, 'b--s','LineWidth',1.8, 'MarkerSize',4);
xlabel('Frazione nodi rimossi f');
ylabel('Giant Component S');
title('Frammentazione della rete');
legend('Rimozione mirata','Rimozione casuale','Location','southwest');
grid on;

subplot(1,2,2);
plot(frac, apl_targ, 'r-o', 'LineWidth',1.8, 'MarkerSize',4); hold on;
plot(frac, apl_rand, 'b--s','LineWidth',1.8, 'MarkerSize',4);
xlabel('Frazione nodi rimossi f');
ylabel('Average path length \langle l \rangle');
title('Efficienza residua della rete');
legend('Rimozione mirata','Rimozione casuale','Location','northwest');
grid on;

% Calcola Giant Component e APL sulla componente connessa più grande
function [gc_size, apl] = network_metrics(G, N_orig)
    bins = conncomp(G);
    counts = histcounts(bins, 1:max(bins)+1);
    [gc_count, gc_label] = max(counts);
    gc_size = gc_count / N_orig;
    Gsub = subgraph(G, bins == gc_label);
    D = distances(Gsub, 'Method','unweighted');
    mask = ~isinf(D) & D > 0;
    apl = mean(D(mask));
end

% Coefficiente di clustering locale per ogni nodo
function C = local_clustering_coefficient(G)
    N = numnodes(G);
    A = adjacency(G);
    C = zeros(N, 1);
    for i = 1:N
        nb = find(A(i,:));
        k  = length(nb);
        if k < 2; continue; end
        C(i) = sum(sum(A(nb, nb))) / 2 / (k*(k-1)/2);
    end
end
