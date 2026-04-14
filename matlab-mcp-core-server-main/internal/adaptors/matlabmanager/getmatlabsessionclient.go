// Copyright 2025 The MathWorks, Inc.

package matlabmanager

import (
	"context"
	"fmt"

	"github.com/matlab/matlab-mcp-core-server/internal/entities"
)

func (m *MATLABManager) GetMATLABSessionClient(ctx context.Context, sessionLogger entities.Logger, sessionID entities.SessionID) (entities.MATLABSessionClient, error) {
	client, err := m.sessionStore.Get(sessionID)
	if err != nil {
		return nil, err
	}

	pingResponse := client.Ping(ctx, sessionLogger)
	if !pingResponse.IsAlive {
		return nil, fmt.Errorf("MATLAB session %v is not alive", sessionID)
	}

	return client, nil
}
