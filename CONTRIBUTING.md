# Contributing to Spacecraft Self-Shadowing Analysis

Thank you for your interest in contributing! This document provides guidelines for contributing to this project.

## 🚀 Getting Started

1. **Fork the repository** on GitHub
2. **Clone your fork** locally:
   ```bash
   git clone https://github.com/yourusername/Comparative-Analysis-Spacecraft-Self-Shadowing.git
   cd Comparative-Analysis-Spacecraft-Self-Shadowing
   ```
3. **Create a branch** for your changes:
   ```bash
   git checkout -b feature/your-feature-name
   ```

## 📝 Types of Contributions

We welcome various types of contributions:

### 🐛 Bug Reports
- Use the GitHub issue tracker
- Describe the bug clearly with reproduction steps
- Include system information (OS, compiler version, etc.)
- Provide minimal code examples if applicable

### 💡 Feature Requests
- Open an issue describing the feature
- Explain the use case and benefits
- Discuss implementation approach if you have ideas

### 🔧 Code Contributions
- Implement new self-shadowing algorithms
- Improve performance of existing algorithms
- Add visualization tools
- Enhance documentation
- Write tests

## 🛠️ Development Guidelines

### Code Style

#### C++
- Follow modern C++17 standards
- Use meaningful variable names
- Add comments for complex algorithms
- Keep functions focused and concise
- Use const correctness

```cpp
// Good example
std::vector<int> calculate_labels_algorithm_name(
    const std::vector<Triangle>& triangles,
    const std::vector<double>& sun_vector
) {
    // Clear implementation with comments
}
```

#### Python
- Follow PEP 8 style guide
- Use type hints where appropriate
- Write docstrings for functions
- Keep code modular and reusable

```python
def visualize_results(data: pd.DataFrame, output_path: str) -> None:
    """
    Visualize benchmark results.
    
    Args:
        data: DataFrame containing benchmark data
        output_path: Path to save visualization
    """
    pass
```

### Testing

- Test new algorithms against the reference (pairwise) implementation
- Ensure algorithms produce correct results on edge cases
- Verify performance improvements don't sacrifice accuracy
- Test on multiple dataset scales

### Benchmarking

When adding new algorithms:
1. Implement the algorithm following the standard interface
2. Add it to the benchmark suite in `SM3D.cpp`
3. Run comprehensive benchmarks across all datasets
4. Document performance characteristics
5. Update README with results

## 📊 Submitting Changes

1. **Commit your changes** with clear messages:
   ```bash
   git commit -m "Add sweep line algorithm optimization"
   ```

2. **Push to your fork**:
   ```bash
   git push origin feature/your-feature-name
   ```

3. **Create a Pull Request**:
   - Provide a clear title and description
   - Reference related issues
   - Include benchmark results for algorithm changes
   - Ensure all tests pass
   - Update documentation as needed

### Pull Request Checklist

- [ ] Code follows project style guidelines
- [ ] Comments and documentation are updated
- [ ] New algorithms include benchmark results
- [ ] No compiler warnings
- [ ] README updated if needed
- [ ] CHANGELOG updated (if applicable)

## 🔬 Algorithm Contributions

When submitting a new self-shadowing algorithm:

### Required Information
1. **Algorithm name and description**
2. **Time complexity** (theoretical)
3. **Space complexity**
4. **Benchmark results** on standard datasets
5. **Use cases** where it excels
6. **Limitations** or edge cases
7. **References** to papers/sources (if applicable)

### Implementation Requirements
- Follow the standard function signature
- Include inline comments for complex logic
- Add brief description at function start
- Handle edge cases gracefully
- Return consistent output format

### Benchmark Template
```cpp
// Algorithm: Your Algorithm Name
// Time Complexity: O(?)
// Space Complexity: O(?)
// Best for: [describe use case]
// Reference: [paper/source if applicable]

std::vector<int> calculate_labels_your_algorithm(
    const std::vector<Triangle>& triangles,
    const std::vector<double>& sun_vector
) {
    // Implementation
}
```

## 📚 Documentation

- Update README.md for major features
- Add inline code comments
- Include usage examples
- Document algorithm complexity
- Provide benchmark comparisons

## 🤝 Code of Conduct

- Be respectful and constructive
- Welcome newcomers
- Focus on technical merit
- Collaborate openly
- Give credit appropriately

## 💬 Questions?

- Open an issue for general questions
- Tag maintainers for urgent matters
- Check existing issues first
- Provide context in your questions

## 📄 License

By contributing, you agree that your contributions will be licensed under the MIT License.

---

Thank you for contributing to spacecraft dynamics research! 🛰️
