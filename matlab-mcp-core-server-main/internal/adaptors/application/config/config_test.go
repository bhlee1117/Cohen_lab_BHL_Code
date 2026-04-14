// Copyright 2025-2026 The MathWorks, Inc.

package config_test

import (
	"encoding/json"
	"path/filepath"
	"testing"
	"time"

	"github.com/matlab/matlab-mcp-core-server/internal/adaptors/application/config"
	"github.com/matlab/matlab-mcp-core-server/internal/adaptors/application/parameter/defaultparameters"
	"github.com/matlab/matlab-mcp-core-server/internal/entities"
	"github.com/matlab/matlab-mcp-core-server/internal/messages"
	"github.com/matlab/matlab-mcp-core-server/internal/testutils"
	configmocks "github.com/matlab/matlab-mcp-core-server/mocks/adaptors/application/config"
	"github.com/stretchr/testify/assert"
	"github.com/stretchr/testify/require"
)

func defaultParameters() []entities.Parameter {
	return []entities.Parameter{
		defaultparameters.HelpMode(),
		defaultparameters.VersionMode(),
		defaultparameters.DisableTelemetry(),
		defaultparameters.UseSingleMATLABSession(),
		defaultparameters.LogLevel(),
		defaultparameters.PreferredLocalMATLABRoot(),
		defaultparameters.PreferredMATLABStartingDirectory(),
		defaultparameters.BaseDir(),
		defaultparameters.WatchdogMode(),
		defaultparameters.ServerInstanceID(),
		defaultparameters.InitializeMATLABOnStartup(),
		defaultparameters.MATLABDisplayMode(),
		defaultparameters.TelemetryCollectorEndpoint(),
		defaultparameters.TelemetryCollectionInterval(),
		defaultparameters.TelemetryCollectorEndpointInsecure(),
	}
}

func configDefaultParsedArgs() map[string]any {
	result := make(map[string]any)
	for _, p := range defaultParameters() {
		result[p.GetID()] = p.GetDefaultValue()
	}
	return result
}

func TestNewConfig_InvalidParameterType(t *testing.T) {
	testCases := []struct {
		name         string
		key          string
		invalidValue any
		expectedType string
	}{
		{name: "LogLevel wrong type", key: defaultparameters.LogLevel().GetID(), invalidValue: 123, expectedType: "string"},
		{name: "UseSingleMATLABSession wrong type", key: defaultparameters.UseSingleMATLABSession().GetID(), invalidValue: "true", expectedType: "bool"},
		{name: "InitializeMATLABOnStartup wrong type", key: defaultparameters.InitializeMATLABOnStartup().GetID(), invalidValue: "false", expectedType: "bool"},
		{name: "VersionMode wrong type", key: defaultparameters.VersionMode().GetID(), invalidValue: "false", expectedType: "bool"},
		{name: "HelpMode wrong type", key: defaultparameters.HelpMode().GetID(), invalidValue: "false", expectedType: "bool"},
		{name: "DisableTelemetry wrong type", key: defaultparameters.DisableTelemetry().GetID(), invalidValue: "false", expectedType: "bool"},
		{name: "PreferredLocalMATLABRoot wrong type", key: defaultparameters.PreferredLocalMATLABRoot().GetID(), invalidValue: 123, expectedType: "string"},
		{name: "PreferredMATLABStartingDirectory wrong type", key: defaultparameters.PreferredMATLABStartingDirectory().GetID(), invalidValue: 123, expectedType: "string"},
		{name: "BaseDir wrong type", key: defaultparameters.BaseDir().GetID(), invalidValue: 123, expectedType: "string"},
		{name: "WatchdogMode wrong type", key: defaultparameters.WatchdogMode().GetID(), invalidValue: "false", expectedType: "bool"},
		{name: "ServerInstanceID wrong type", key: defaultparameters.ServerInstanceID().GetID(), invalidValue: 123, expectedType: "string"},
		{name: "MATLABDisplayMode wrong type", key: defaultparameters.MATLABDisplayMode().GetID(), invalidValue: 123, expectedType: "string"},
		{name: "TelemetryCollectorEndpoint wrong type", key: defaultparameters.TelemetryCollectorEndpoint().GetID(), invalidValue: 123, expectedType: "string"},
		{name: "TelemetryCollectionInterval wrong type", key: defaultparameters.TelemetryCollectionInterval().GetID(), invalidValue: "1m", expectedType: "time.Duration"},
		{name: "TelemetryCollectorEndpointInsecure wrong type", key: defaultparameters.TelemetryCollectorEndpointInsecure().GetID(), invalidValue: "false", expectedType: "bool"},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			// Arrange
			mockOSLayer := &configmocks.MockOSLayer{}
			defer mockOSLayer.AssertExpectations(t)

			mockParser := &configmocks.MockParser{}
			defer mockParser.AssertExpectations(t)

			mockBuildInfo := &configmocks.MockBuildInfo{}
			defer mockBuildInfo.AssertExpectations(t)

			programName := "testprocess"
			args := []string{programName}

			parsedArgs := configDefaultParsedArgs()
			parsedArgs[tc.key] = tc.invalidValue

			mockOSLayer.EXPECT().
				Args().
				Return(args).
				Once()

			mockParser.EXPECT().
				Parse(args[1:]).
				Return([]entities.Parameter{}, parsedArgs, []string{}, nil).
				Once()

			expectedError := messages.New_StartupErrors_InvalidParameterType_Error(tc.key, tc.expectedType)

			// Act
			cfg, err := config.NewConfig(mockOSLayer, mockParser, mockBuildInfo)

			// Assert
			require.Equal(t, expectedError, err)
			assert.Nil(t, cfg)
		})
	}
}

func TestNewConfig_MissingParameter(t *testing.T) {
	testCases := []struct {
		name       string
		missingKey string
	}{
		{name: "missing LogLevel", missingKey: defaultparameters.LogLevel().GetID()},
		{name: "missing UseSingleMATLABSession", missingKey: defaultparameters.UseSingleMATLABSession().GetID()},
		{name: "missing InitializeMATLABOnStartup", missingKey: defaultparameters.InitializeMATLABOnStartup().GetID()},
		{name: "missing VersionMode", missingKey: defaultparameters.VersionMode().GetID()},
		{name: "missing HelpMode", missingKey: defaultparameters.HelpMode().GetID()},
		{name: "missing DisableTelemetry", missingKey: defaultparameters.DisableTelemetry().GetID()},
		{name: "missing PreferredLocalMATLABRoot", missingKey: defaultparameters.PreferredLocalMATLABRoot().GetID()},
		{name: "missing PreferredMATLABStartingDirectory", missingKey: defaultparameters.PreferredMATLABStartingDirectory().GetID()},
		{name: "missing BaseDir", missingKey: defaultparameters.BaseDir().GetID()},
		{name: "missing WatchdogMode", missingKey: defaultparameters.WatchdogMode().GetID()},
		{name: "missing ServerInstanceID", missingKey: defaultparameters.ServerInstanceID().GetID()},
		{name: "missing MATLABDisplayMode", missingKey: defaultparameters.MATLABDisplayMode().GetID()},
		{name: "missing TelemetryCollectorEndpoint", missingKey: defaultparameters.TelemetryCollectorEndpoint().GetID()},
		{name: "missing TelemetryCollectionInterval", missingKey: defaultparameters.TelemetryCollectionInterval().GetID()},
		{name: "missing TelemetryCollectorEndpointInsecure", missingKey: defaultparameters.TelemetryCollectorEndpointInsecure().GetID()},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			// Arrange
			mockOSLayer := &configmocks.MockOSLayer{}
			defer mockOSLayer.AssertExpectations(t)

			mockParser := &configmocks.MockParser{}
			defer mockParser.AssertExpectations(t)

			mockBuildInfo := &configmocks.MockBuildInfo{}
			defer mockBuildInfo.AssertExpectations(t)

			programName := "testprocess"
			args := []string{programName}

			parsedArgs := configDefaultParsedArgs()
			delete(parsedArgs, tc.missingKey)

			mockOSLayer.EXPECT().
				Args().
				Return(args).
				Once()

			mockParser.EXPECT().
				Parse(args[1:]).
				Return([]entities.Parameter{}, parsedArgs, []string{}, nil).
				Once()

			expectedError := messages.New_StartupErrors_InvalidParameterKey_Error(tc.missingKey)

			// Act
			cfg, err := config.NewConfig(mockOSLayer, mockParser, mockBuildInfo)

			// Assert
			require.Equal(t, expectedError, err)
			assert.Nil(t, cfg)
		})
	}
}

func TestNewConfig_ParseError(t *testing.T) {
	// Arrange
	mockOSLayer := &configmocks.MockOSLayer{}
	defer mockOSLayer.AssertExpectations(t)

	mockParser := &configmocks.MockParser{}
	defer mockParser.AssertExpectations(t)

	mockBuildInfo := &configmocks.MockBuildInfo{}
	defer mockBuildInfo.AssertExpectations(t)

	programName := "testprocess"
	args := []string{programName}

	mockOSLayer.EXPECT().
		Args().
		Return(args).
		Once()

	mockParser.EXPECT().
		Parse(args[1:]).
		Return(nil, nil, nil, messages.AnError).
		Once()

	// Act
	cfg, err := config.NewConfig(mockOSLayer, mockParser, mockBuildInfo)

	// Assert
	require.ErrorIs(t, err, messages.AnError)
	assert.Nil(t, cfg, "Config should be nil")
}

func TestNewConfig_InvalidLogLevel(t *testing.T) {
	// Arrange
	mockOSLayer := &configmocks.MockOSLayer{}
	defer mockOSLayer.AssertExpectations(t)

	mockParser := &configmocks.MockParser{}
	defer mockParser.AssertExpectations(t)

	mockBuildInfo := &configmocks.MockBuildInfo{}
	defer mockBuildInfo.AssertExpectations(t)

	programName := "testprocess"
	args := []string{programName}

	parsedArgs := configDefaultParsedArgs()
	parsedArgs[defaultparameters.LogLevel().GetID()] = "invalid-level"

	mockOSLayer.EXPECT().
		Args().
		Return(args).
		Once()

	mockParser.EXPECT().
		Parse(args[1:]).
		Return([]entities.Parameter{}, parsedArgs, []string{}, nil).
		Once()

	expectedError := messages.New_StartupErrors_InvalidLogLevel_Error("invalid-level")

	// Act
	cfg, err := config.NewConfig(mockOSLayer, mockParser, mockBuildInfo)

	// Assert
	require.Equal(t, expectedError, err)
	assert.Nil(t, cfg, "Config should be nil")
}

func TestConfig_Version_HappyPath(t *testing.T) {
	// Arrange
	mockOSLayer := &configmocks.MockOSLayer{}
	defer mockOSLayer.AssertExpectations(t)

	mockParser := &configmocks.MockParser{}
	defer mockParser.AssertExpectations(t)

	mockBuildInfo := &configmocks.MockBuildInfo{}
	defer mockBuildInfo.AssertExpectations(t)

	programName := "testprocess"
	args := []string{programName}

	expectedVersion := "github.com/matlab/matlab-mcp-core-server v1.2.3"

	mockOSLayer.EXPECT().
		Args().
		Return(args).
		Once()

	mockParser.EXPECT().
		Parse(args[1:]).
		Return([]entities.Parameter{}, configDefaultParsedArgs(), []string{}, nil).
		Once()

	mockBuildInfo.EXPECT().
		FullVersion().
		Return(expectedVersion).
		Once()

	// Act
	cfg, err := config.NewConfig(mockOSLayer, mockParser, mockBuildInfo)
	require.NoError(t, err)

	version := cfg.Version()

	// Assert
	require.Equal(t, expectedVersion, version)
}

func TestConfig_SpecifiedParameters_HappyPath(t *testing.T) {
	// Arrange
	mockOSLayer := &configmocks.MockOSLayer{}
	defer mockOSLayer.AssertExpectations(t)

	mockParser := &configmocks.MockParser{}
	defer mockParser.AssertExpectations(t)

	mockBuildInfo := &configmocks.MockBuildInfo{}
	defer mockBuildInfo.AssertExpectations(t)

	programName := "testprocess"
	args := []string{programName}
	expectedSpecifiedParameters := []string{"DisableTelemetry", "LogLevel"}

	mockOSLayer.EXPECT().
		Args().
		Return(args).
		Once()

	mockParser.EXPECT().
		Parse(args[1:]).
		Return(defaultParameters(), configDefaultParsedArgs(), expectedSpecifiedParameters, nil).
		Once()

	// Act
	cfg, err := config.NewConfig(mockOSLayer, mockParser, mockBuildInfo)
	require.NoError(t, err)

	result := cfg.SpecifiedParameters()

	// Assert
	require.Equal(t, expectedSpecifiedParameters, result)
}

func TestConfig_InitializeMATLABOnStartup_DisabledWhenNotSingleSession(t *testing.T) {
	// Arrange
	mockOSLayer := &configmocks.MockOSLayer{}
	defer mockOSLayer.AssertExpectations(t)

	mockParser := &configmocks.MockParser{}
	defer mockParser.AssertExpectations(t)

	mockBuildInfo := &configmocks.MockBuildInfo{}
	defer mockBuildInfo.AssertExpectations(t)

	programName := "testprocess"
	args := []string{programName}

	parsedArgs := configDefaultParsedArgs()
	parsedArgs[defaultparameters.UseSingleMATLABSession().GetID()] = false
	parsedArgs[defaultparameters.InitializeMATLABOnStartup().GetID()] = true

	mockOSLayer.EXPECT().
		Args().
		Return(args).
		Once()

	mockParser.EXPECT().
		Parse(args[1:]).
		Return([]entities.Parameter{}, parsedArgs, []string{}, nil).
		Once()

	// Act
	cfg, err := config.NewConfig(mockOSLayer, mockParser, mockBuildInfo)

	// Assert
	require.NoError(t, err)
	assert.False(t, cfg.InitializeMATLABOnStartup(), "InitializeMATLABOnStartup should be false when UseSingleMATLABSession is false")
}

func TestNewConfig_TelemetryCollectionInterval_FallsBackToDefaultWhenNotPositive(t *testing.T) {
	testCases := []struct {
		name     string
		interval time.Duration
	}{
		{name: "zero interval", interval: 0},
		{name: "negative interval", interval: -time.Minute},
	}

	for _, tc := range testCases {
		t.Run(tc.name, func(t *testing.T) {
			// Arrange
			mockOSLayer := &configmocks.MockOSLayer{}
			defer mockOSLayer.AssertExpectations(t)

			mockParser := &configmocks.MockParser{}
			defer mockParser.AssertExpectations(t)

			mockBuildInfo := &configmocks.MockBuildInfo{}
			defer mockBuildInfo.AssertExpectations(t)

			programName := "testprocess"
			args := []string{programName}

			parsedArgs := configDefaultParsedArgs()
			parsedArgs[defaultparameters.TelemetryCollectionInterval().GetID()] = tc.interval

			expectedInterval := defaultparameters.TelemetryCollectionInterval().GetTypedDefaultValue()

			mockOSLayer.EXPECT().
				Args().
				Return(args).
				Once()

			mockParser.EXPECT().
				Parse(args[1:]).
				Return([]entities.Parameter{}, parsedArgs, []string{}, nil).
				Once()

			// Act
			cfg, err := config.NewConfig(mockOSLayer, mockParser, mockBuildInfo)

			// Assert
			require.NoError(t, err)
			assert.Equal(t, expectedInterval, cfg.TelemetryCollectionInterval())
		})
	}
}

func TestConfig_RecordToLogger_HappyPath(t *testing.T) {
	// Arrange
	parsedArgs := configDefaultParsedArgs()
	parsedArgs[defaultparameters.DisableTelemetry().GetID()] = true
	parsedArgs[defaultparameters.PreferredMATLABStartingDirectory().GetID()] = filepath.Join("home", "user")
	parsedArgs[defaultparameters.LogLevel().GetID()] = string(entities.LogLevelDebug)
	parsedArgs[defaultparameters.PreferredLocalMATLABRoot().GetID()] = filepath.Join("home", "matlab")
	parsedArgs[defaultparameters.UseSingleMATLABSession().GetID()] = false

	expectedLogMessage := "Configuration state"
	expectedConfigField := map[string]any{
		defaultparameters.DisableTelemetry().GetID():                 true,
		defaultparameters.PreferredMATLABStartingDirectory().GetID(): filepath.Join("home", "user"),
		defaultparameters.LogLevel().GetID():                         string(entities.LogLevelDebug),
		defaultparameters.PreferredLocalMATLABRoot().GetID():         filepath.Join("home", "matlab"),
		defaultparameters.UseSingleMATLABSession().GetID():           false,
		defaultparameters.InitializeMATLABOnStartup().GetID():        false,
	}

	parameters := []entities.Parameter{
		defaultparameters.DisableTelemetry(),
		defaultparameters.UseSingleMATLABSession(),
		defaultparameters.LogLevel(),
		defaultparameters.PreferredLocalMATLABRoot(),
		defaultparameters.PreferredMATLABStartingDirectory(),
		defaultparameters.InitializeMATLABOnStartup(),
		defaultparameters.HelpMode(),
		defaultparameters.VersionMode(),
		defaultparameters.BaseDir(),
		defaultparameters.WatchdogMode(),
		defaultparameters.ServerInstanceID(),
	}

	mockOSLayer := &configmocks.MockOSLayer{}
	defer mockOSLayer.AssertExpectations(t)

	mockParser := &configmocks.MockParser{}
	defer mockParser.AssertExpectations(t)

	mockBuildInfo := &configmocks.MockBuildInfo{}
	defer mockBuildInfo.AssertExpectations(t)

	programName := "testprocess"
	args := []string{programName}

	mockParser.EXPECT().
		Parse(args[1:]).
		Return(parameters, parsedArgs, []string{}, nil)

	mockOSLayer.EXPECT().
		Args().
		Return(args).
		Once()

	cfg, err := config.NewConfig(mockOSLayer, mockParser, mockBuildInfo)
	require.NoError(t, err)

	testLogger := testutils.NewInspectableLogger()

	// Act
	cfg.RecordToLogger(testLogger)

	// Assert
	infoLogs := testLogger.InfoLogs()
	require.Len(t, infoLogs, 1)

	fields, found := infoLogs[expectedLogMessage]
	require.True(t, found, "Expected log message not found")

	for expectedField, expectedValue := range expectedConfigField {
		actualValue, exists := fields[expectedField]
		require.True(t, exists, "%s field not found in log", expectedField)
		assert.Equal(t, expectedValue, actualValue, "%s field has incorrect value", expectedField)
	}
}

func TestConfig_AsPIISafeJSONString_HappyPath(t *testing.T) {
	// Arrange
	mockOSLayer := &configmocks.MockOSLayer{}
	defer mockOSLayer.AssertExpectations(t)

	mockParser := &configmocks.MockParser{}
	defer mockParser.AssertExpectations(t)

	mockBuildInfo := &configmocks.MockBuildInfo{}
	defer mockBuildInfo.AssertExpectations(t)

	programName := "testprocess"
	args := []string{programName}

	parameters := defaultParameters()
	parsedArgs := configDefaultParsedArgs()

	mockOSLayer.EXPECT().
		Args().
		Return(args).
		Once()

	mockParser.EXPECT().
		Parse(args[1:]).
		Return(parameters, parsedArgs, []string{}, nil).
		Once()

	// Act
	cfg, err := config.NewConfig(mockOSLayer, mockParser, mockBuildInfo)
	require.NoError(t, err)

	result := cfg.AsPIISafeJSONString()

	// Assert
	var parsed map[string]any
	require.NoError(t, json.Unmarshal([]byte(result), &parsed))

	piiSafeParams := []entities.Parameter{
		defaultparameters.HelpMode(),
		defaultparameters.VersionMode(),
		defaultparameters.DisableTelemetry(),
		defaultparameters.UseSingleMATLABSession(),
		defaultparameters.LogLevel(),
		defaultparameters.InitializeMATLABOnStartup(),
		defaultparameters.MATLABDisplayMode(),
		defaultparameters.WatchdogMode(),
		defaultparameters.TelemetryCollectionInterval(),
		defaultparameters.TelemetryCollectorEndpointInsecure(),
	}
	for _, param := range piiSafeParams {
		var expected any
		raw, _ := json.Marshal(parsedArgs[param.GetID()])
		_ = json.Unmarshal(raw, &expected)
		assert.Equal(t, expected, parsed[param.GetID()], "%s should show actual value", param.GetID())
	}

	redactedParams := []entities.Parameter{
		defaultparameters.PreferredLocalMATLABRoot(),
		defaultparameters.PreferredMATLABStartingDirectory(),
		defaultparameters.BaseDir(),
		defaultparameters.ServerInstanceID(),
		defaultparameters.TelemetryCollectorEndpoint(),
	}
	for _, param := range redactedParams {
		assert.Equal(t, config.RedactedValue, parsed[param.GetID()], "%s should be redacted", param.GetID())
	}
}
