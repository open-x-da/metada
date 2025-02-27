// Unified Traits Class that combines all dependencies
template <typename... Backends>
struct AppTraits {};

// Specialization for Logger and Config backends
template <typename LoggerBackend, typename ConfigBackend>
struct AppTraits<LoggerBackend, ConfigBackend> {
    using LoggerType = LoggerBackend;            // Logger backend
    using ConfigType = ConfigBackend;            // Config backend
};

// Future specializations can be added as needed:
// template <typename LoggerBackend, typename ConfigBackend, typename NewBackend>
// struct AppTraits<LoggerBackend, ConfigBackend, NewBackend> { ... };
