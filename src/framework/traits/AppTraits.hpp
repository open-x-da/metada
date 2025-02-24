// Unified Traits Class that combines all dependencies
template <typename LoggerBackend, typename ConfigBackend>
struct AppTraits {
    using LoggerType = LoggerBackend;            // Logger backend
    using ConfigType = ConfigBackend;            // Config backend
};
