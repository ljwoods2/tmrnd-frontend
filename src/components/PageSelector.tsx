interface PageSelectorProps {
  currentPage: number;
  onPageChange: (page: number) => void;
}

const PageSelector: React.FC<PageSelectorProps> = ({
  currentPage,
  onPageChange,
}) => (
  <div className="page-selector">
    {[0, 1, 2].map((page) => (
      <span
        key={page}
        className={`dot ${currentPage === page ? "active" : ""}`}
        onClick={() => onPageChange(page)}
      />
    ))}
  </div>
);

export default PageSelector;
