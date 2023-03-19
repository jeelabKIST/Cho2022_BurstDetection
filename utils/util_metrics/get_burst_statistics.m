function [bstat_class] = get_burst_statistics(btime_true, btrue_cycl, btime_test, dt, foi)
    %% Function 'get_burst_statistics'
    % DESCRIPTION
    % Extracts classification statistics from the burst simulations and
    % resulting detections
    
    % INPUT
    %    Variable       Data Type           Description
    % 1. btime_true     [1 x N cell]      : set of time arrays of simluated bursts
    % 2. btrue_cycl     [1 x N cell]      : set of burst lengths of simulated bursts
    % 3. btime_test     [1 x N cell]      : set of time arrays of detected bursts
    % 4. dt             [rational Q > 0]  : time step by which the window slides (seconds)
    % 5. foi            [integer Z > 0]   : frequency of interest based on which the bursts were detected
    
    % OUTPUT
    %    Variable       Data Type           Description
    % 1. bstat_class    [N x N cell]      : table consisting of burst IDs of simulated and detected bursts and classification labels of each detection
    
    % Written by SungJun Cho, August 15, 2021
    % Last Modified on February 25, 2023
    %% Preallocate Output of Arbitrarily Large Size
    truth_class = cell(1,length(btime_true)*length(btime_test));
    %% True Positives
    k = 1;
    for refIdx = 1:length(btime_true)
        alpha = btime_true{refIdx}(1);
        omega = btime_true{refIdx}(end);
        for testIdx = 1:length(btime_test)
            begin = btime_test{testIdx}(1);
            conclude = btime_test{testIdx}(end);
            % Pass if there is no overlap between true and detected bursts
            if begin > omega || conclude < alpha
                continue;
            end
            % Case: Detected burst conatins true burst OR complete overlap
            if begin <= alpha && conclude >= omega
                refID = refIdx; bID = testIdx;
                bClass = 'tp';
                cycle = btrue_cycl(refIdx);
                ta = compute_temporal_accuracy(alpha,omega,begin,conclude);
                truth_class{k} = {refID; bID; bClass; cycle; ta};
                k = k + 1;
            else
                % Case: True burst contains detected burst
                if begin > alpha && conclude < omega
                    refID = refIdx; bID = testIdx;
                    bClass = 'tp';
                    cycle = btrue_cycl(refIdx);
                    ta = compute_temporal_accuracy(alpha,omega,begin,conclude);
                    truth_class{k} = {refID; bID; bClass; cycle; ta};
                    k = k + 1;
                % Case: Partial overlap in which detected burst precedes true burst
                elseif begin <= alpha && conclude < omega
                    front_ratio = (conclude - alpha) / (conclude - begin);
                    if front_ratio >= 0.3 && front_ratio <= 1
                        refID = refIdx; bID = testIdx;
                        bClass = 'tp';
                        cycle = btrue_cycl(refIdx);
                        ta = compute_temporal_accuracy(alpha,omega,begin,conclude);
                        truth_class{k} = {refID; bID; bClass; cycle; ta};
                        k = k + 1;
                    end
                % Case: Partial overlap in which detected burst succeeds true burst
                elseif conclude >= omega && begin > alpha
                    back_ratio = (omega - begin) / (conclude - begin);
                    if back_ratio >= 0.3 && back_ratio <= 1
                        refID = refIdx; bID = testIdx;
                        bClass = 'tp';
                        cycle = btrue_cycl(refIdx);
                        ta = compute_temporal_accuracy(alpha,omega,begin,conclude);
                        truth_class{k} = {refID; bID; bClass; cycle; ta};
                        k = k + 1;
                    end
                end
            end
        end
    end
    % Confirm preallocation was large enough
    if k > size(truth_class,2)
        warning(['The size of current output exceeded the size preallocated at the outset. ' ...
            'This can happen when no bursts were detected by the algorithm.']);
    end
    % Reorganize array containing statistics of true positives
    truth_class = [truth_class{:}];
    %% False Positives
    % Case: Detected bursts with no overlapping true bursts
    if isempty(truth_class)
        set_bID = 0;                                    % null set of true positive bursts (if none exists)
    else
        set_bID = unique(cell2mat(truth_class(2,:)));   % numerical set of true positive bursts
    end
    set_bID_all = 1:length(btime_test);                 % numerical set of all detected bursts
    fpIdx = ~ismember(set_bID_all,set_bID);             % get indices of falsely detected bursts
    fpList = set_bID_all(fpIdx);                        % numerical set of false positive bursts
    if ~isempty(fpList)
        num_cycles = cellfun(@(x) num2cell(round(length(x)*(foi*dt))), btime_test(fpList));
        false_positive_class = [num2cell(NaN(1,length(fpList))); num2cell(fpList);
                                repmat({'fp'},1,length(fpList)); num_cycles;
                                num2cell(NaN(1,length(fpList)))];
        % Output: Class of True and False Positive Bursts
        bstat_class = cat(2,truth_class,false_positive_class);
    else
        bstat_class = truth_class;
    end
    %% False Negatives
    % Case: True bursts with no overlapping detected bursts
    if isempty(truth_class)
        set_refID = 0;                                    % null set of true reference bursts (if none exists)
    else
        set_refID = unique(cell2mat(truth_class(1,:)));   % numerical set of true reference bursts
    end
    set_refID_all = 1:length(btime_true);                 % nuermical set of all reference bursts
    fnIdx = ~ismember(set_refID_all,set_refID);           % get indices of falsely detected bursts
    fnList = set_refID_all(fnIdx);                        % numerical set of flase negative bursts
    if ~isempty(fnList)
        false_negative_class = [num2cell(fnList); num2cell(NaN(1,length(fnList)));
                                repmat({'fn'},1,length(fnList)); num2cell(btrue_cycl(fnIdx));
                                num2cell(NaN(1,length(fnList)))];
        % Output: Class of True Positive, False Positive, and False Negative Bursts
        bstat_class = cat(2,bstat_class,false_negative_class);
    end
end

%% Appendix: In-Script Functions
% Function #1: Computes temporal concurrence of a single burst
function [ta] = compute_temporal_accuracy(s_start,s_stop,d_start,d_stop)
    numerator = min(s_stop,d_stop) - max(s_start,d_start);
    denominator = max(s_stop,d_stop) - min(s_start,d_start);
    ta = numerator / denominator;
end