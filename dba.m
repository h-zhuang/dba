% sequence_set: a cell array with sequences to be averaged
sequence_set = {[],[]};
% num_iter: number of iterations
num_iter = 1;
avg = dtw_avg(sequence_set,num_iter);

function avg = dtw_avg(sequence_set,num_iter)
    
    num_seq = length(sequence_set);
    s = sequence_set{1};
    max_length = length(s);
    for j = 1:1:num_seq
        if length(sequence_set{j}) >= max_length
            max_length = length(sequence_set{j});
            longest_seq = sequence_set{j};
            max_length_ind = j;
        end
    end
    average_init = longest_seq;
    for i = 1:1:num_iter
        avg = dba_update(average_init,sequence_set,max_length_ind);
    end
end

function updated_average = dba_update(average_init,sequence_set,max_length_ind)
    num_seq = length(sequence_set);
    updated_average = average_init;
    for i = 1:1:num_seq
        if i ~= max_length_ind
            s = sequence_set{i};
            s_alignment = dtw_alignment(updated_average,s);
            updated_average = (updated_average + s_alignment)/2;
        end
    end
end


function s_aligned = dtw_alignment(sref,s)
    cost = cost_matrix(sref,s);
    l = length(sref);
    s_aligned = zeros([1,l]);
    i = size(cost,1);
    j = size(cost,2);
    while(i>=1) && (j>=1)
        s_aligned(i) = s(j);
        if i == 1 
            j = j-1;
        elseif j == 1 
            i = i-1;
        else
            score = min(cost(i-1,j-1),min(cost(i,j-1),cost(i-1,j)));
            if score == cost(i-1,j-1)
                i = i-1;
                j = j-1;
            elseif score == cost(i-1,j)
                i = i-1;
            else
                j = j-1;
            end
        end
    end
end

function cost = cost_matrix(sref,s)
% Generate the cost matrix
cost = zeros([length(sref),length(s)]);
    for i = 1:1:length(sref)
        for j = 1:1:length(s)
            cost(i,j) = abs(sref(i)-s(j));
        end
    end
end