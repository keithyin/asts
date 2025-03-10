build:
	cargo build --release

install:
	cp target/release/asts /usr/bin/
	cp target/release/asrtc /usr/bin/

clean:
	rm -rf target