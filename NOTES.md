# Spike Result: diffsol compilation on wasm32-unknown-unknown

## Summary
The attempt to compile `diffsol` as part of `crates/core` for the `wasm32-unknown-unknown` target failed. 

## Technical Details
- **Dependency**: `diffsol` (v0.14) transitively depends on `rand` / `getrandom` (v0.3.4).
- **Error**: `getrandom v0.3.4` does not support `wasm32-unknown-unknown` by default:
  ```
  error: The wasm32-unknown-unknown targets are not supported by default; you may need to enable the "wasm_js" configuration flag.
  ```
- **Architectural Constraints**: Since `crates/core` is a plain library crate with **no wasm dependencies** (an architectural boundary to prevent any imports of `wasm-bindgen` or `web-sys`), we cannot easily configure custom JS-backed random generators or target-specific features in `core`.

## Decision / Action Taken
Following the task instructions:
1. Removed `diffsol` from the dependencies of `crates/core`.
2. Added `diffsol` exclusively to the dependencies of `crates/reference-cli` (which runs natively and has no wasm32 constraints).
