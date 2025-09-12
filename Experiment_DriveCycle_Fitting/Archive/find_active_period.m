function [start_idx, end_idx] = find_active_period(I_data, current_threshold)
    % 유효한 전류 구간 찾기
    % Input:
    %   I_data: 전류 데이터 [A]
    %   current_threshold: 전류 임계값 [A] (기본값: 0.1A)
    % Output:
    %   start_idx: 유효 구간 시작 인덱스
    %   end_idx: 유효 구간 끝 인덱스
    
    if nargin < 2
        current_threshold = 0.1;  % 기본 임계값 0.1A
    end
    
    % 절대값이 임계값 이상인 구간 찾기
    active_indices = find(abs(I_data) >= current_threshold);
    
    if isempty(active_indices)
        % 유효한 구간이 없으면 전체 구간 사용
        start_idx = 1;
        end_idx = length(I_data);
        fprintf('Warning: No active period found. Using entire dataset.\n');
    else
        % 첫 번째와 마지막 유효 인덱스
        start_idx = active_indices(1);
        end_idx = active_indices(end);
        
        % 최소 1000개 포인트는 확보
        if (end_idx - start_idx + 1) < 1000
            center = round((start_idx + end_idx) / 2);
            start_idx = max(1, center - 500);
            end_idx = min(length(I_data), center + 500);
        end
        
        fprintf('Active period: %d ~ %d (%.1f%% of total)\n', ...
            start_idx, end_idx, (end_idx - start_idx + 1) / length(I_data) * 100);
    end
end 