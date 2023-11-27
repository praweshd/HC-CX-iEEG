function [matchedref,matched2,unmatchedref,unmatched2,match_idxref,match_idx2] = LFPmatching(res1,res2,tolerance)
%USAGE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%USAGE  
%LFP matching made to identify and isolate IED events in NG (res2) that  % 
%matched/unmatched with events in HP (res1).                             %
% Very similar to IED matching function                                  %
%Input: res1: must be the reference res CH (e.g. Hippocampus)            %
%       res2: must be the CH whose events you want to isolate (spindle)  %
%       tol : tolerance (ms)to match/unmatch two events, default: 20 ms  % 
%                                                                        %
%Output: matched: events in res2 that matched with res1                  %
%        unmatched: events in res2 that didn't match with res1           %
%        match_idx: indices of res2 that matched with res1               %
%                                                                        %
% Created April 8,2020 - Prawesh Dahal (TNL,Columbia)                    %
%------------------------------------------------------------------------%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tol = tolerance/1000;
matched2 = []; unmatched2 = []; match_idx2 = [];
matchedref = []; unmatchedref = []; match_idxref = [];

for i = 1:length(res1)
    
    A=res1(i);
    flag = 0;
    for j = 1:length(res2)
        
        B=res2(j);
        
        if abs(A-B) <= tol
            matched2 = [matched2 B];
            match_idx2 = [match_idx2 j];
            
            matchedref = [matchedref A];
            match_idxref = [match_idxref i];
            
%             disp(['HP=', num2str(i),' NG=', num2str(j)]);
            flag = 1;
            break
        end 
        
        if flag ==1 
            break
        end 
    end 
    
end 

unmatched2 = setdiff(res2,matched2);
unmatchedref = setdiff(res1,matchedref);

matched2 = matched2';
matchedref = matchedref';

end

