build:
	cargo build --release

install:
	cp target/release/asts /usr/bin/

clean:
	rm -rf target