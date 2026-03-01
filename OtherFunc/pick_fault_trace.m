%% pick_fault_trace.m
% Interactive fault trace picker for InSAR .grd data
%
% Instructions:
%   - Click pairs of points on the plot to define fault segments
%   - Each pair of clicks = one fault segment (lon1 lat1 lon2 lat2)
%   - Press ANY KEY to finish picking and optionally save the output
%
% Output format (one segment per line):
%   lon1 lat1 lon2 lat2
%   lon3 lat3 lon4 lat4
%   ...

clear; clc; close all;

%% Update your path here!
addpath(genpath('/Volumes/T7/Research/PamirProject'))

%% --- Load GRD file ---
[fname, fpath] = uigetfile('*.grd', 'Select InSAR .grd file');
if isequal(fname, 0)
    disp('No file selected. Exiting.');
    return;
end
filepath = fullfile(fpath, fname);

fprintf('Loading: %s\n', filepath);
[x, y, z] = grdread2(filepath);

%% --- Plot the InSAR data ---
fig = figure('Name', 'InSAR Fault Trace Picker', 'NumberTitle', 'off');
imagesc(x, y, z);
axis xy;
colormap(jet);
colorbar;
xlabel('Longitude');
ylabel('Latitude');
title({'InSAR Displacement Map — Fault Trace Picker', ...
       'Click pairs of points to define fault segments | Press any key to finish'});
set(gca, 'FontSize', 11);

%% --- Interactive picking loop ---
segments   = [];   % will grow as [lon1 lat1 lon2 lat2; ...]
pt_buffer  = [];   % holds the first click of a pair
seg_lines  = [];   % handles to plotted segment lines
pt_markers = [];   % handles to plotted click markers

fprintf('\n--- Picking mode active ---\n');
fprintf('Click PAIRS of points to define fault segments.\n');
fprintf('Press any key when done.\n\n');

% Use a flag controlled by KeyPressFcn
done = false;
set(fig, 'KeyPressFcn', @(src,evt) setappdata(src, 'done', true));
setappdata(fig, 'done', false);

while ishandle(fig)
    % Check if user pressed a key
    if getappdata(fig, 'done')
        done = true;
        break;
    end

    % Wait for a click (ginput with timeout-style using waitforbuttonpress)
    try
        % waitforbuttonpress returns 0 for mouse click, 1 for key press
        val = waitforbuttonpress;
    catch
        % Figure was closed
        break;
    end

    if ~ishandle(fig)
        break;
    end

    if val == 1
        % Key was pressed — finish
        done = true;
        break;
    end

    % It was a mouse click — get coordinates
    cp = get(gca, 'CurrentPoint');
    cx = cp(1, 1);
    cy = cp(1, 2);

    if isempty(pt_buffer)
        % First click of a pair — buffer it
        pt_buffer = [cx, cy];
        hold on;
        h = plot(cx, cy, 'wo', 'MarkerFaceColor', 'w', ...
                 'MarkerSize', 6, 'LineWidth', 1.2);
        pt_markers = [pt_markers, h];
        fprintf('  Point 1 of segment: (%.5f, %.5f) — click a second point\n', cx, cy);
    else
        % Second click — complete the segment
        x1 = pt_buffer(1);  y1 = pt_buffer(2);
        x2 = cx;            y2 = cy;

        segments = [segments; x1 y1 x2 y2];

        % Plot the segment
        hold on;
        hl = plot([x1 x2], [y1 y2], 'w-', 'LineWidth', 2);
        h2 = plot(x2, y2, 'wo', 'MarkerFaceColor', 'w', ...
                  'MarkerSize', 6, 'LineWidth', 1.2);
        % Label segment number
        seg_num = size(segments, 1);
        text((x1+x2)/2, (y1+y2)/2, sprintf(' %d', seg_num), ...
             'Color', 'w', 'FontSize', 9, 'FontWeight', 'bold');

        seg_lines  = [seg_lines,  hl];
        pt_markers = [pt_markers, h2];

        fprintf('  Segment %d: (%.5f, %.5f) -> (%.5f, %.5f)\n', ...
                seg_num, x1, y1, x2, y2);

        pt_buffer = [];   % reset buffer for next segment
    end
end

%% --- Close figure ---
if ishandle(fig)
    close(fig);
end

%% --- Report results ---
n_seg = size(segments, 1);
fprintf('\n--- Picking complete: %d segment(s) defined ---\n', n_seg);

if n_seg == 0
    fprintf('No segments were picked. Exiting.\n');
    return;
end

% Print summary
fprintf('\nSegments (lon1 lat1 lon2 lat2):\n');
for i = 1:n_seg
    fprintf('  %.6f  %.6f  %.6f  %.6f\n', segments(i,1), segments(i,2), ...
                                           segments(i,3), segments(i,4));
end

%% --- Ask user whether to save ---
choice = questdlg('Save fault trace to .txt file?', ...
                  'Save Output', 'Yes', 'No', 'Yes');

if strcmp(choice, 'Yes')

    % --- Ask for output format ---
    fmt_choice = questdlg( ...
        ['Choose output format:' newline newline ...
         'Segment rows:  lon1 lat1 lon2 lat2' newline ...
         '               lon3 lat3 lon4 lat4' newline newline ...
         'Point rows:    lon1 lat1' newline ...
         '               lon2 lat2' newline ...
         '               lon3 lat3' newline ...
         '               lon4 lat4'], ...
        'Output Format', ...
        'Segment rows (endpoint pairs)', ...
        'Point rows (one point per line)', ...
        'Segment rows (endpoint pairs)');   % default

    if isempty(fmt_choice)
        fprintf('Save cancelled.\n');
        return;
    end

    [out_fname, out_path] = uiputfile('*.txt', 'Save fault trace as', ...
                                      'fault_trace.txt');
    if isequal(out_fname, 0)
        fprintf('Save cancelled.\n');
    else
        out_file = fullfile(out_path, out_fname);
        fid = fopen(out_file, 'w');
        if fid == -1
            error('Could not open file for writing: %s', out_file);
        end

        if strcmp(fmt_choice, 'Segment rows (endpoint pairs)')
            % Format: lon1 lat1 lon2 lat2  (one segment per line)
            for i = 1:n_seg
                fprintf(fid, '%.6f  %.6f  %.6f  %.6f\n', ...
                        segments(i,1), segments(i,2), ...
                        segments(i,3), segments(i,4));
            end
            fprintf('Format: segment rows (lon1 lat1 lon2 lat2)\n');

        else
            % Format: lon lat  (one point per line, two lines per segment)
            for i = 1:n_seg
                fprintf(fid, '%.6f  %.6f\n', segments(i,1), segments(i,2));
                fprintf(fid, '%.6f  %.6f\n', segments(i,3), segments(i,4));
            end
            fprintf('Format: point rows (lon lat, one per line)\n');
        end

        fclose(fid);
        fprintf('Fault trace saved to: %s\n', out_file);
    end
else
    fprintf('Output not saved.\n');
end