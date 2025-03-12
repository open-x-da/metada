// Unified Traits Class that combines all dependencies
template <typename... Backends>
struct AppTraits {};

// Specialization for Logger and Config backends
template <typename LoggerBackend, typename ConfigBackend>
struct AppTraits<LoggerBackend, ConfigBackend> {
    using LoggerType = LoggerBackend;            // Logger backend
    using ConfigType = ConfigBackend;            // Config backend
};

// Specialization for Logger, Config and State backends
template <typename LoggerBackend, typename ConfigBackend, typename StateBackend>
struct AppTraits<LoggerBackend, ConfigBackend, StateBackend> {
    using LoggerType = LoggerBackend;            // Logger backend
    using ConfigType = ConfigBackend;            // Config backend
    using StateType = StateBackend;              // State backend
};

// Specialization for Logger, Config, State and Observation backends
template <typename LoggerBackend, typename ConfigBackend, typename StateBackend, typename ObservationBackend>
struct AppTraits<LoggerBackend, ConfigBackend, StateBackend, ObservationBackend> {
    using LoggerType = LoggerBackend;            // Logger backend
    using ConfigType = ConfigBackend;            // Config backend
    using StateType = StateBackend;              // State backend
    using ObservationType = ObservationBackend; // Observation backend
};

// Specialization for Logger, Config, State, Observation and ObsOperator backends
template <typename LoggerBackend, typename ConfigBackend, typename StateBackend, typename ObservationBackend, typename ObsOperatorBackend>
struct AppTraits<LoggerBackend, ConfigBackend, StateBackend, ObservationBackend, ObsOperatorBackend> {
    using LoggerType = LoggerBackend;            // Logger backend
    using ConfigType = ConfigBackend;            // Config backend
    using StateType = StateBackend;              // State backend
    using ObservationType = ObservationBackend; // Observation backend
    using ObsOperatorType = ObsOperatorBackend; // ObsOperator backend
};

// Specialization for Logger, Config, State, Observation, ObsOperator and Model backends
template <typename LoggerBackend, typename ConfigBackend, typename StateBackend, typename ObservationBackend, typename ObsOperatorBackend, typename ModelBackend>
struct AppTraits<LoggerBackend, ConfigBackend, StateBackend, ObservationBackend, ObsOperatorBackend, ModelBackend> {
    using LoggerType = LoggerBackend;            // Logger backend
    using ConfigType = ConfigBackend;            // Config backend
    using StateType = StateBackend;              // State backend
    using ObservationType = ObservationBackend; // Observation backend
    using ObsOperatorType = ObsOperatorBackend; // ObsOperator backend
    using ModelType = ModelBackend;               // Model backend
};

// Specialization for Logger, Config, State, Observation, ObsOperator and Model backends
template <typename LoggerBackend, typename ConfigBackend, typename GeometryBackend, typename StateBackend, typename ObservationBackend, typename ObsOperatorBackend, typename ModelBackend>
struct AppTraits<LoggerBackend, ConfigBackend, GeometryBackend, StateBackend, ObservationBackend, ObsOperatorBackend, ModelBackend> {
    using LoggerType = LoggerBackend;            // Logger backend
    using ConfigType = ConfigBackend;            // Config backend
    using GeometryType = GeometryBackend;        // Geometry backend
    using StateType = StateBackend;              // State backend
    using ObservationType = ObservationBackend; // Observation backend
    using ObsOperatorType = ObsOperatorBackend; // ObsOperator backend
    using ModelType = ModelBackend;               // Model backend
};