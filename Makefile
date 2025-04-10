build:
	cargo build --release

install:
	cargo build --release
	cp target/release/asts /usr/bin/
	cp target/release/asrtc /usr/bin/
	cp target/release/asfstc /usr/bin/

clean:
	rm -rf target