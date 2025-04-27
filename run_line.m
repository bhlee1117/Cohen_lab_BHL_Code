function run_line(filename, sectionIndex)
    % Read entire script
    lines = readlines(filename, TextType="string");

    % Find section delimiters (lines that start with %%)
    sectionStartIdx = [1; find(startsWith(strtrim(lines), "%%")) + 1];
    sectionEndIdx = [sectionStartIdx(2:end) - 1; numel(lines)];

    % Validate input
    if sectionIndex < 1 || sectionIndex > numel(sectionStartIdx)
        error('Invalid section index. Script has %d sections.', numel(sectionStartIdx));
    end

    % Extract and run the selected section
    codeBlock = lines(sectionStartIdx(sectionIndex):sectionEndIdx(sectionIndex));
    codeStr = join(codeBlock, newline);
    evalin('base', codeStr);
end