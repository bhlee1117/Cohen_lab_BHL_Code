// Copyright 2025 The MathWorks, Inc.

package responseconverter

import (
	"github.com/matlab/matlab-mcp-core-server/internal/adaptors/mcp/tools"
	"github.com/matlab/matlab-mcp-core-server/internal/entities"
)

func ConvertEvalResponseToRichContent(response entities.EvalResponse) tools.RichContent {
	imageData := make([]tools.PNGImageData, len(response.Images))
	for i := range response.Images {
		imageData[i] = tools.PNGImageData(response.Images[i])
	}
	return tools.RichContent{
		TextContent:  []string{response.ConsoleOutput},
		ImageContent: imageData,
	}
}
