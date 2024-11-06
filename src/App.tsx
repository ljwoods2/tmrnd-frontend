import { useState } from "react";
import "./App.css";
import SubmitJobs from "./components/SubmitJobs.tsx";
import PageSelector from "./components/PageSelector.tsx";
import TestQueue from "./components/TestQueue.tsx";
import AccessInstructions from "./components/AccessInstructions.tsx";

const App: React.FC = () => {
  const [currentPage, setCurrentPage] = useState(0);

  const renderPage = () => {
    switch (currentPage) {
      case 0:
        return <SubmitJobs />;
      case 1:
        return <TestQueue />;
      case 2:
        return <AccessInstructions />;
      default:
        return <h1>Error</h1>;
    }
  };

  return (
    <div className="background-gradient">
      <div className="content">
        <div className="container">{renderPage()}</div>
        <PageSelector currentPage={currentPage} onPageChange={setCurrentPage} />
      </div>
    </div>
  );
};

export default App;
