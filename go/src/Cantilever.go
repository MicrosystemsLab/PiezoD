package main
import "fmt"
import "strconv"

type Cantilever struct {
	cantileverType int
}

func NewCantilever() *Cantilever {
	return &Cantilever{cantileverType: 2}
}

func PrintCantilever(c Cantilever) (string) {
	return strconv.Itoa(c.cantileverType)
}

func main() {
	var c Cantilever = *NewCantilever()

	fmt.Printf("Hello world!\n")
	fmt.Printf(PrintCantilever(c))
}
