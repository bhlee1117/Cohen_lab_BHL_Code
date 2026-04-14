// Copyright 2026 The MathWorks, Inc.

package selector_test

import (
	"testing"

	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/mock"

	"github.com/matlab/matlab-mcp-core-server/internal/adaptors/application/definition"
	"github.com/matlab/matlab-mcp-core-server/internal/adaptors/application/parameter/defaultparameters/selector"
	"github.com/matlab/matlab-mcp-core-server/internal/messages"
	selectormocks "github.com/matlab/matlab-mcp-core-server/mocks/adaptors/application/parameter/defaultparameters/selector"
)

func TestSelector_DefaultParameters_DescriptionsResolved(t *testing.T) {
	// Arrange
	mockAppDef := &selectormocks.MockApplicationDefinition{}
	defer mockAppDef.AssertExpectations(t)

	mockMessageCatalog := &selectormocks.MockMessageCatalog{}
	defer mockMessageCatalog.AssertExpectations(t)

	expectedDescriptions := map[messages.MessageKey]string{
		messages.CLIMessages_HelpDescription:                             "Help description",
		messages.CLIMessages_VersionDescription:                          "Version description",
		messages.CLIMessages_DisableTelemetryDescription:                 "Disable telemetry description",
		messages.CLIMessages_BaseDirDescription:                          "Base dir description",
		messages.CLIMessages_LogLevelDescription:                         "Log level description",
		messages.CLIMessages_InternalUseDescription:                      "Internal use description",
		messages.CLIMessages_PreferredLocalMATLABRootDescription:         "MATLAB root description",
		messages.CLIMessages_PreferredMATLABStartingDirectoryDescription: "MATLAB starting directory description",
		messages.CLIMessages_UseSingleMATLABSessionDescription:           "Single MATLAB session description",
		messages.CLIMessages_InitializeMATLABOnStartupDescription:        "Initialize MATLAB on startup description",
		messages.CLIMessages_DisplayModeDescription:                      "Display mode description",
	}

	mockAppDef.EXPECT().
		Features().
		Return(definition.Features{}).
		Once()

	for key, description := range expectedDescriptions {
		times := 1
		if key == messages.CLIMessages_InternalUseDescription {
			times = 5
		}

		mockMessageCatalog.EXPECT().
			Get(key).
			Return(description).
			Times(times)
	}

	sut := selector.New(mockAppDef, mockMessageCatalog)

	// Act
	parameters := sut.DefaultParameters()

	// Assert
	for _, p := range parameters {
		assert.NotEmpty(t, p.GetDescription(), "parameter %s should have a description", p.GetID())
	}
}

func TestSelector_DefaultParameters_MATLABEnabled(t *testing.T) {
	// Arrange
	mockAppDef := &selectormocks.MockApplicationDefinition{}
	defer mockAppDef.AssertExpectations(t)

	mockMessageCatalog := &selectormocks.MockMessageCatalog{}
	defer mockMessageCatalog.AssertExpectations(t)

	mockAppDef.EXPECT().
		Features().
		Return(definition.Features{MATLAB: definition.MATLABFeature{Enabled: true}}).
		Once()

	mockMessageCatalog.EXPECT().
		Get(mock.Anything).
		Return("description")

	sut := selector.New(mockAppDef, mockMessageCatalog)

	// Act
	parameters := sut.DefaultParameters()

	// Assert
	assert.Len(t, parameters, 15)

	for _, p := range parameters {
		assert.True(t, p.GetActive(), "parameter %s should be active", p.GetID())
	}
}

func TestSelector_DefaultParameters_MATLABDisabled(t *testing.T) {
	// Arrange
	mockAppDef := &selectormocks.MockApplicationDefinition{}
	defer mockAppDef.AssertExpectations(t)

	mockMessageCatalog := &selectormocks.MockMessageCatalog{}
	defer mockMessageCatalog.AssertExpectations(t)

	commonParameterIDs := map[string]bool{
		"HelpMode":                           true,
		"VersionMode":                        true,
		"DisableTelemetry":                   true,
		"BaseDir":                            true,
		"LogLevel":                           true,
		"WatchdogMode":                       true,
		"ServerInstanceID":                   true,
		"TelemetryCollectorEndpoint":         true,
		"TelemetryCollectionInterval":        true,
		"TelemetryCollectorEndpointInsecure": true,
	}

	matlabParameterIDs := map[string]bool{
		"PreferredLocalMATLABRoot":         true,
		"PreferredMATLABStartingDirectory": true,
		"UseSingleMATLABSession":           true,
		"InitializeMATLABOnStartup":        true,
		"MATLABDisplayMode":                true,
	}

	mockAppDef.EXPECT().
		Features().
		Return(definition.Features{MATLAB: definition.MATLABFeature{Enabled: false}}).
		Once()

	mockMessageCatalog.EXPECT().
		Get(mock.Anything).
		Return("description")

	sut := selector.New(mockAppDef, mockMessageCatalog)

	// Act
	parameters := sut.DefaultParameters()

	// Assert
	assert.Len(t, parameters, 15)

	for _, p := range parameters {
		if commonParameterIDs[p.GetID()] {
			assert.True(t, p.GetActive(), "common parameter %s should be active", p.GetID())
		}
		if matlabParameterIDs[p.GetID()] {
			assert.False(t, p.GetActive(), "MATLAB parameter %s should be inactive", p.GetID())
		}
	}
}
