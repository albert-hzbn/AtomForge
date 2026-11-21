//! Linux: the library finds spglib beside it (the portable package layout).
fn main() {
    if std::env::var("TARGET").unwrap_or_default().contains("linux") {
        println!("cargo:rustc-link-arg-cdylib=-Wl,-rpath,$ORIGIN");
    }
}
